"""Storage backend abstraction for pipeline run results.

Supports two backends:
    - ``LocalStorage``: reads/writes files on the local filesystem under
      ``results/`` (development mode).
    - ``SupabaseStorage``: reads/writes files to a Supabase Storage bucket
      (production mode).

The factory ``get_storage()`` returns the appropriate backend based on
environment variables.  Upload operations are always best-effort — errors
are logged but never crash the pipeline.

Usage:
    from core.storage import get_storage

    storage = get_storage()
    storage.upload_file("pipeline-results", "run123/output.csv", Path("out.csv"))
    data = storage.download_file("pipeline-results", "run123/output.csv")
"""

from __future__ import annotations

import mimetypes
import os
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

from core.logger import get_logger

logger = get_logger("storage")

# Default bucket name, configurable via env var
STORAGE_BUCKET: str = os.getenv("STORAGE_BUCKET", "pipeline-results")


# ---------------------------------------------------------------------------
# Abstract base
# ---------------------------------------------------------------------------


class StorageBackend(ABC):
    """Abstract interface for file storage operations."""

    @abstractmethod
    def upload_file(self, bucket: str, path: str, file_path: Path) -> str:
        """Upload a local file to storage.

        Args:
            bucket: Storage bucket / namespace.
            path: Destination path within the bucket.
            file_path: Local file to upload.

        Returns:
            The storage path of the uploaded file.
        """

    @abstractmethod
    def upload_bytes(
        self, bucket: str, path: str, data: bytes, content_type: str
    ) -> str:
        """Upload raw bytes to storage.

        Args:
            bucket: Storage bucket / namespace.
            path: Destination path within the bucket.
            data: Raw byte content.
            content_type: MIME type (e.g. ``application/json``).

        Returns:
            The storage path of the uploaded file.
        """

    @abstractmethod
    def download_file(self, bucket: str, path: str) -> bytes:
        """Download a file's contents from storage.

        Args:
            bucket: Storage bucket / namespace.
            path: Path within the bucket.

        Returns:
            The file contents as bytes.

        Raises:
            FileNotFoundError: If the file does not exist.
        """

    @abstractmethod
    def list_files(self, bucket: str, prefix: str) -> list[str]:
        """List files under a given prefix.

        Args:
            bucket: Storage bucket / namespace.
            prefix: Path prefix to filter by.

        Returns:
            List of file paths matching the prefix.
        """

    @abstractmethod
    def get_public_url(self, bucket: str, path: str) -> str:
        """Return a publicly-accessible URL for a stored file.

        Args:
            bucket: Storage bucket / namespace.
            path: Path within the bucket.

        Returns:
            A URL string.
        """

    def upload_directory(
        self, bucket: str, prefix: str, local_dir: Path
    ) -> list[str]:
        """Upload all files in a directory tree to storage.

        Files are uploaded with paths relative to ``local_dir`` appended
        to ``prefix``.  Errors on individual files are logged and skipped.

        Args:
            bucket: Storage bucket / namespace.
            prefix: Base path prefix in the bucket.
            local_dir: Local directory to upload.

        Returns:
            List of storage paths that were successfully uploaded.
        """
        uploaded: list[str] = []
        if not local_dir.exists() or not local_dir.is_dir():
            logger.warning("upload_directory: %s does not exist or is not a directory", local_dir)
            return uploaded

        for file_path in sorted(local_dir.rglob("*")):
            if not file_path.is_file():
                continue
            relative = file_path.relative_to(local_dir)
            dest_path = f"{prefix}/{relative}".replace("\\", "/")
            try:
                self.upload_file(bucket, dest_path, file_path)
                uploaded.append(dest_path)
            except Exception:
                logger.warning(
                    "Failed to upload %s to %s/%s (skipping)",
                    file_path,
                    bucket,
                    dest_path,
                )
        return uploaded


# ---------------------------------------------------------------------------
# Local filesystem backend (development)
# ---------------------------------------------------------------------------


class LocalStorage(StorageBackend):
    """Stores files on the local filesystem under a base directory.

    In local mode the ``bucket`` parameter is used as a subdirectory
    name inside ``base_dir``.  This mirrors the bucket/path structure
    of remote storage so paths are interchangeable.
    """

    def __init__(self, base_dir: Path | None = None) -> None:
        if base_dir is None:
            base_dir = Path(__file__).resolve().parent.parent / "results"
        self.base_dir = base_dir

    def _resolve(self, bucket: str, path: str) -> Path:
        """Resolve a bucket/path pair to a local filesystem path."""
        resolved = (self.base_dir / bucket / path).resolve()
        # Prevent path traversal
        base_resolved = self.base_dir.resolve()
        if not str(resolved).startswith(str(base_resolved)):
            raise ValueError("Path traversal detected")
        return resolved

    def upload_file(self, bucket: str, path: str, file_path: Path) -> str:
        dest = self._resolve(bucket, path)
        dest.parent.mkdir(parents=True, exist_ok=True)

        # Copy file contents
        dest.write_bytes(file_path.read_bytes())
        logger.debug("LocalStorage: copied %s -> %s", file_path, dest)
        return f"{bucket}/{path}"

    def upload_bytes(
        self, bucket: str, path: str, data: bytes, content_type: str
    ) -> str:
        dest = self._resolve(bucket, path)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(data)
        logger.debug("LocalStorage: wrote %d bytes to %s", len(data), dest)
        return f"{bucket}/{path}"

    def download_file(self, bucket: str, path: str) -> bytes:
        resolved = self._resolve(bucket, path)
        if not resolved.exists():
            raise FileNotFoundError(f"File not found: {bucket}/{path}")
        return resolved.read_bytes()

    def list_files(self, bucket: str, prefix: str) -> list[str]:
        base = self._resolve(bucket, prefix)
        if not base.exists():
            return []
        results: list[str] = []
        if base.is_file():
            return [f"{bucket}/{prefix}"]
        for file_path in sorted(base.rglob("*")):
            if file_path.is_file():
                relative = file_path.relative_to(self.base_dir / bucket)
                results.append(str(relative).replace("\\", "/"))
        return results

    def get_public_url(self, bucket: str, path: str) -> str:
        # Local storage has no public URL — return a file:// URI
        resolved = self._resolve(bucket, path)
        return resolved.as_uri()


# ---------------------------------------------------------------------------
# Supabase Storage backend (production)
# ---------------------------------------------------------------------------


class SupabaseStorage(StorageBackend):
    """Stores files in a Supabase Storage bucket.

    Uses the lazily-initialised Supabase client from ``core.db``.
    Path format inside the bucket: ``{team_id}/{run_id}/{filename}``
    or ``{run_id}/{filename}`` when no team_id is available.
    """

    def __init__(self) -> None:
        self._client: Any | None = None

    @property
    def client(self) -> Any:
        """Lazy-load the Supabase client."""
        if self._client is None:
            from core.db import _get_client
            self._client = _get_client()
        return self._client

    def upload_file(self, bucket: str, path: str, file_path: Path) -> str:
        content_type = mimetypes.guess_type(str(file_path))[0] or "application/octet-stream"
        data = file_path.read_bytes()
        return self.upload_bytes(bucket, path, data, content_type)

    def upload_bytes(
        self, bucket: str, path: str, data: bytes, content_type: str
    ) -> str:
        self.client.storage.from_(bucket).upload(
            path=path,
            file=data,
            file_options={"content-type": content_type, "upsert": "true"},
        )
        logger.info("SupabaseStorage: uploaded %s/%s (%d bytes)", bucket, path, len(data))
        return f"{bucket}/{path}"

    def download_file(self, bucket: str, path: str) -> bytes:
        response = self.client.storage.from_(bucket).download(path)
        if response is None:
            raise FileNotFoundError(f"File not found in Supabase: {bucket}/{path}")
        return response

    def list_files(self, bucket: str, prefix: str) -> list[str]:
        response = self.client.storage.from_(bucket).list(path=prefix)
        if not response:
            return []
        results: list[str] = []
        for item in response:
            name = item.get("name", "")
            if name:
                full_path = f"{prefix}/{name}".lstrip("/")
                results.append(full_path)
        return results

    def get_public_url(self, bucket: str, path: str) -> str:
        result = self.client.storage.from_(bucket).get_public_url(path)
        return result


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

_storage_instance: StorageBackend | None = None


def get_storage() -> StorageBackend:
    """Return the appropriate storage backend.

    Uses ``SupabaseStorage`` when ``SUPABASE_URL`` is set in the
    environment, otherwise falls back to ``LocalStorage``.

    The instance is cached as a module-level singleton.
    """
    global _storage_instance  # noqa: PLW0603
    if _storage_instance is not None:
        return _storage_instance

    supabase_url = os.getenv("SUPABASE_URL", "")
    if supabase_url:
        try:
            _storage_instance = SupabaseStorage()
            logger.info("Using SupabaseStorage backend (bucket=%s)", STORAGE_BUCKET)
        except Exception:
            logger.warning(
                "Failed to initialise SupabaseStorage, falling back to LocalStorage"
            )
            _storage_instance = LocalStorage()
    else:
        _storage_instance = LocalStorage()
        logger.info("Using LocalStorage backend")

    return _storage_instance


def build_storage_prefix(run_id: str, team_id: str | None = None) -> str:
    """Build the storage path prefix for a run's outputs.

    Returns ``{team_id}/{run_id}`` when a team_id is available,
    otherwise ``{run_id}``.
    """
    if team_id:
        return f"{team_id}/{run_id}"
    return run_id
