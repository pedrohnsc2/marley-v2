"""RQ worker entry point for pipeline execution.

Usage (inside Docker container):
    python -m core.worker

This starts an rq worker that listens on the "pipelines" queue
and executes pipeline runs in the background. Each job runs
_execute_run_job() from worker_api.py.

The worker process is long-lived and handles one job at a time.
Scale by running more replicas in docker-compose.
"""

from __future__ import annotations

import os
import logging

from dotenv import load_dotenv

load_dotenv()

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s %(message)s",
)
logger = logging.getLogger("worker")


def main() -> None:
    import redis
    from rq import Queue, Worker

    redis_url = os.getenv("REDIS_URL", "redis://localhost:6379/0")
    conn = redis.from_url(redis_url)

    # Verify Redis is reachable
    try:
        conn.ping()
        logger.info("Connected to Redis at %s", redis_url)
    except redis.ConnectionError:
        logger.error("Cannot connect to Redis at %s", redis_url)
        raise

    queue = Queue("pipelines", connection=conn)
    worker = Worker([queue], connection=conn, name=f"marley-worker-{os.getpid()}")

    logger.info("Starting Marley pipeline worker (PID %d)", os.getpid())
    logger.info("Listening on queue: pipelines")

    worker.work(with_scheduler=False)


if __name__ == "__main__":
    main()
