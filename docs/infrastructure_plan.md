# Marley — Infrastructure Plan for Multi-User Scale

> Consolidado por 4 agentes: Architect, Tech Lead, PO, Infra (2026-04-14)

---

## Executive Summary

Scale Marley from 1 researcher to 10-50 across UFMG, Fiocruz, and CTVacinas.
Total estimated effort: **24-34 dev days** across 5 phases.
Production host: **Hostinger KVM 1** (already provisioned, active until 2027-05).

---

## Architecture Target State

```
Internet
  |
  v
[Caddy] :443 automatic HTTPS (Let's Encrypt)
  |
  +-- /* ---------> [Next.js web] :3000
  +-- /api/runs/start -> [Worker API (FastAPI)] :8000
  |                          |
  |                     [Redis] :6379
  |                          |
  |               [Worker 1] [Worker 2]  (rq pipeline workers)
  |                          |
  +-- Supabase Realtime <----+
  |                          |
  v                          v
[Supabase]            [Supabase Storage]
  - Auth                - pipeline-results/
  - PostgreSQL            {team_id}/{run_id}/
  - RLS on all tables
```

---

## Phase A: Authentication Foundation (3-4 days)

### Database Migrations

```sql
-- 002_user_profiles.sql
create table profiles (
  id uuid primary key references auth.users(id) on delete cascade,
  email text not null,
  full_name text not null default '',
  institution text default '',
  orcid text,
  role text not null default 'researcher'
    check (role in ('admin', 'pi', 'researcher', 'viewer')),
  created_at timestamptz not null default now(),
  updated_at timestamptz not null default now()
);

-- Auto-create profile on signup
create function handle_new_user() returns trigger as $$
begin
  insert into profiles (id, email, full_name)
  values (new.id, new.email, coalesce(new.raw_user_meta_data->>'full_name', ''));
  return new;
end;
$$ language plpgsql security definer;

create trigger on_auth_user_created
  after insert on auth.users
  for each row execute function handle_new_user();

-- 003_runs_add_user_id.sql
alter table pipeline_runs add column user_id uuid references auth.users(id);
create index idx_pipeline_runs_user_id on pipeline_runs(user_id);
```

### Files to Create
- `web/lib/supabase/client.ts` — browser client with `createBrowserClient`
- `web/lib/supabase/server.ts` — server client with cookie handling
- `web/lib/supabase/middleware.ts` — session refresh
- `web/app/[locale]/login/page.tsx` — email/password + magic link
- `web/app/[locale]/signup/page.tsx` — registration
- `web/app/auth/callback/route.ts` — OAuth callback

### Files to Modify
- `web/middleware.ts` — add JWT validation, protect routes
- `web/lib/api-auth.ts` — replace static key with Supabase session
- `web/app/api/runs/start/route.ts` — extract user_id from session
- `core/run.py` — add user_id to PipelineRun
- `core/launcher.py` — pass user_id
- `core/db_runs.py` — include user_id in upsert

### Auth Methods (priority order)
1. Email + magic link (lowest friction)
2. Email + password (fallback)
3. ORCID OAuth (researchers already have IDs)

---

## Phase B: Authorization & Multi-tenancy (5-7 days)

### Database Migrations

```sql
-- 004_teams.sql
create table teams (
  id uuid primary key default gen_random_uuid(),
  name text not null,
  slug text unique not null,
  institution text default '',
  settings jsonb not null default '{}',
  created_at timestamptz not null default now()
);

create table team_memberships (
  team_id uuid not null references teams(id) on delete cascade,
  user_id uuid not null references auth.users(id) on delete cascade,
  role text not null default 'researcher'
    check (role in ('owner', 'admin', 'researcher', 'viewer')),
  accepted_at timestamptz,
  primary key (team_id, user_id)
);

alter table pipeline_runs add column team_id uuid references teams(id);

-- 005_rls_all_tables.sql
-- RLS helper functions
create function user_team_ids() returns uuid[] as $$
  select array_agg(team_id) from team_memberships
  where user_id = auth.uid() and accepted_at is not null
$$ language sql security definer stable;

-- Apply RLS to pipeline_runs, pipeline_stages, and all 6 science tables
-- Pattern: team_id = any(user_team_ids()) or user_id is null (legacy)

-- 006_api_keys.sql
create table api_keys (
  id uuid primary key default gen_random_uuid(),
  user_id uuid not null references auth.users(id),
  team_id uuid not null references teams(id),
  key_hash text not null,
  key_prefix char(8) not null unique,
  name text not null,
  scopes text[] not null default '{run:create,run:read}',
  expires_at timestamptz,
  created_at timestamptz not null default now()
);
```

### Role Permissions

| Action | owner | admin | researcher | viewer |
|--------|-------|-------|-----------|--------|
| Create runs | yes | yes | yes | no |
| View team runs | yes | yes | yes | yes |
| Manage members | yes | yes | no | no |
| Create presets | yes | yes | yes | no |
| API key mgmt | yes | yes | own only | no |

---

## Phase C: Scalable Pipeline Execution (7-10 days)

### Job Queue: rq + Redis

Replace `execFileSync` + `subprocess.Popen` with proper job queue.

**Why rq over Celery**: simpler, same Redis, zero config, adequate for 2-5 concurrent runs.
**Why not pgmq**: avoids Postgres extension dependency, Redis is $0-5/month.

```python
# core/worker.py
from rq import Worker, Queue
from redis import Redis

redis_conn = Redis.from_url(os.environ['REDIS_URL'])
queue = Queue('pipeline-jobs', connection=redis_conn)

def process_pipeline_job(run_id: str):
    mgr = RunManager()
    run = mgr.load_run(run_id)
    execute_pipeline_run(run, run.parameters)

if __name__ == '__main__':
    worker = Worker([queue], connection=redis_conn)
    worker.work()
```

### Worker API (FastAPI)

```python
# core/worker_api.py
@app.post("/api/runs/start")
async def start_run(body: RunRequest):
    run = create_pipeline_run(body.pipeline, body.parameters)
    queue.enqueue(process_pipeline_job, run.run_id, job_timeout=1800)
    return {"run_id": run.run_id, "status": "queued"}
```

### Storage Migration

Results move from local filesystem to Supabase Storage:
- Worker uploads outputs after run completion
- Frontend reads from Supabase Storage (with local fallback for dev)
- StorageBackend abstraction in `web/lib/storage-loader.ts`

---

## Phase D: Collaborative Features (5-7 days)

### User Presets

```sql
-- 009_presets.sql
create table presets (
  id uuid primary key default gen_random_uuid(),
  pipeline text not null,
  name text not null,
  scope text not null default 'user' check (scope in ('system','team','user')),
  team_id uuid references teams(id),
  created_by uuid references auth.users(id),
  display_name text not null,
  description text default '',
  recommended boolean default false,
  parameters jsonb not null,
  forked_from uuid references presets(id),
  created_at timestamptz not null default now()
);
```

### Resolution Order
1. User-scoped preset (this user + team)
2. Team-scoped preset (shared by team)
3. System preset (the 18 JSON files, imported to DB)

### Run Sharing
- Generate signed token per run (30-day validity)
- Read-only link: `/shared/{token}` renders results without auth
- Fork button: clone config to own account

### Publication Export
- "Copy Methods Text" button — generates paragraph for papers
- BibTeX entry with run_id + git_sha
- Reproducibility ZIP: params + results + metadata

---

## Phase E: Deployment & Operations (4-6 days)

### Hosting: Hostinger KVM 1 (already provisioned)

- **OS**: Ubuntu 24.04 LTS
- **CPU**: 1 vCPU (x86_64)
- **RAM**: 4 GB (48% in use at idle, ~2 GB available)
- **Disk**: 50 GB NVMe SSD (14 GB used, 36 GB available)
- **Bandwidth**: 4 TB/month
- **SSH**: `root@<server-ip>`
- **Expires**: 2027-05-27
- **Firewall**: 0 rules configured (needs setup)
- **Malware detector**: Active

#### Capacity Constraints (KVM 1)

With 2 GB free RAM, the stack must be lean:

| Service | RAM Budget | CPU | Notes |
|---------|-----------|-----|-------|
| Caddy | 30 MB | minimal | Reverse proxy + auto HTTPS |
| Next.js | 300 MB | 0.3 | Standalone build (~150 MB image) |
| Redis | 50 MB | minimal | maxmemory 50mb, allkeys-lru |
| Worker (1x) | 1-2 GB | 1.0 | Single concurrent pipeline |
| **Total** | **~1.5-2.5 GB** | ~1.3 | Fits within 2 GB headroom |

**Limitations**: Only 1 pipeline worker (not 2). Runs execute sequentially if queue has multiple jobs. Adequate for 10-20 users with typical usage patterns.

**Upgrade path**: Hostinger KVM 2 (2 vCPU, 8 GB RAM, 100 GB) allows 2-3 concurrent workers for heavier load.

### Docker Compose (production — tuned for KVM 1)

```yaml
services:
  caddy:       # Reverse proxy + auto HTTPS
    mem_limit: 64m
  web:         # Next.js standalone
    mem_limit: 512m
    cpus: 0.5
  worker-api:  # FastAPI (co-located with web to save RAM)
    mem_limit: 256m
    cpus: 0.3
  worker:      # rq pipeline worker (1 replica only on KVM 1)
    mem_limit: 2g
    cpus: 1.0
    deploy:
      replicas: 1
  redis:       # Job queue
    mem_limit: 64m
    command: redis-server --maxmemory 50mb --maxmemory-policy allkeys-lru
```

#### Firewall Rules (to configure on Hostinger)

| Port | Protocol | Source | Service |
|------|----------|--------|---------|
| 22 | TCP | your IP only | SSH |
| 80 | TCP | any | HTTP (Caddy redirects to 443) |
| 443 | TCP | any | HTTPS (Caddy) |
| All others | — | blocked | — |

### CI/CD: GitHub Actions
- `ci.yml`: pytest + next build + docker build on every PR
- `deploy.yml`: push images to ghcr.io, SSH pull on VPS, docker compose up

### Monitoring
- Tier 1 (free): Uptime Robot + Sentry free tier
- Tier 2 (optional): Prometheus + Grafana + cAdvisor in Docker

### Cost Estimate (with Hostinger KVM 1)

| Component | Cost |
|-----------|------|
| Hostinger KVM 1 (already paid through 2027-05) | ~$5/mo (pre-paid) |
| Supabase free tier (500 MB DB, 1 GB storage) | $0 |
| Domain (if custom) | ~$3/mo |
| GitHub Actions + GHCR (public repo) | $0 |
| Sentry free tier | $0 |
| **Total** | **~$5-8/mo** |

**Scaling**: If >20 concurrent users, upgrade to KVM 2 (~$10/mo) for 2x CPU/RAM.

---

## Implementation Timeline

| Phase | Depends on | Days | Key Deliverable |
|-------|-----------|------|-----------------|
| A: Auth | — | 3-4 | Users can login, runs owned by user |
| B: AuthZ | A | 5-7 | Teams, roles, RLS, API keys |
| C: Execution | B | 7-10 | Job queue, remote workers, storage |
| D: Collaboration | B | 5-7 | User presets, sharing, export |
| E: Deployment | C+D | 4-6 | Docker, CI/CD, monitoring |
| **Total** | | **24-34** | **Multi-user production** |

---

## Product MoSCoW (v1.0)

### MUST
1. Auth/accounts (Supabase Auth)
2. Run history with search/filter
3. Preset management (save from UI)
4. Result export (CSV + PDF)
5. Reproducibility metadata

### SHOULD
6. Team/lab model
7. Run sharing (read-only links)
8. Notifications (run complete/failed)
9. Preset sharing/forking

### COULD
10. Comments on results
11. Publication workflow (BibTeX, CITATION.cff)
12. Audit trail
13. Role-based access

### WON'T (v1)
- Multi-tenancy isolation
- Custom pipeline builder
- Billing/quotas
- SSO/SAML
- Mobile app

---

## Success Metrics

| Metric | Target (3 months) |
|--------|-------------------|
| Active users/week | 10+ |
| Runs/week | 30+ |
| Time to first run | < 10 min |
| Run success rate | > 85% |
| Institutions | 3+ |

**North Star**: Runs cited in published papers.

---

## ADR Summary

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Multi-tenancy | Shared schema + RLS | 10-50 teams, schema-per-tenant is overkill |
| Job queue | rq + Redis | Simple, cheap, adequate for scale |
| Auth | Supabase Auth | Free, JWT+RLS integration, email+magic link |
| Deployment | Docker Compose on VPS | K8s premature for <100 users |
| Storage | Supabase Storage | S3-compatible, included in plan, same auth |
| Reverse proxy | Caddy | Auto HTTPS, zero maintenance |
| Hosting | Hostinger KVM 1 (already provisioned, Ubuntu 24.04) | Pre-paid through 2027, adequate for 10-20 users |
