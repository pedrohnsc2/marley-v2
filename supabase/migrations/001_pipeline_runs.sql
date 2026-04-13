-- Pipeline Runs: tracks each execution of a bioinformatics pipeline
-- Created for Marley dynamic pipeline system

-- Table: pipeline_runs
create table if not exists pipeline_runs (
  run_id text primary key,
  pipeline text not null,
  status text not null default 'created'
    check (status in ('created', 'running', 'completed', 'failed', 'cancelled')),
  created_at timestamptz not null default now(),
  started_at timestamptz,
  completed_at timestamptz,
  git_sha text default '',
  parameters jsonb not null default '{}',
  total_duration_s double precision not null default 0,
  output_dir text default '',
  notes text default '',
  tags text[] default '{}',
  updated_at timestamptz not null default now()
);

-- Auto-update updated_at on row change
create or replace function update_updated_at()
returns trigger as $$
begin
  new.updated_at = now();
  return new;
end;
$$ language plpgsql;

create trigger pipeline_runs_updated_at
  before update on pipeline_runs
  for each row execute function update_updated_at();

-- Indexes for common queries
create index if not exists idx_pipeline_runs_pipeline on pipeline_runs(pipeline);
create index if not exists idx_pipeline_runs_status on pipeline_runs(status);
create index if not exists idx_pipeline_runs_created_at on pipeline_runs(created_at desc);

-- Table: pipeline_stages
create table if not exists pipeline_stages (
  id bigint generated always as identity primary key,
  run_id text not null references pipeline_runs(run_id) on delete cascade,
  stage_id text not null,
  name text not null,
  status text not null default 'pending'
    check (status in ('pending', 'running', 'success', 'failed', 'skipped')),
  started_at timestamptz,
  completed_at timestamptz,
  duration_s double precision not null default 0,
  error text,
  key_metrics jsonb not null default '{}',
  constraint unique_run_stage unique (run_id, stage_id)
);

-- Index for stage lookups by run
create index if not exists idx_pipeline_stages_run_id on pipeline_stages(run_id);

-- Row Level Security
alter table pipeline_runs enable row level security;
alter table pipeline_stages enable row level security;

-- Anon users can read (for the web frontend with NEXT_PUBLIC key)
create policy "anon_read_runs" on pipeline_runs
  for select using (true);

create policy "anon_read_stages" on pipeline_stages
  for select using (true);

-- Only service_role can write (Python backend uses service key)
create policy "service_write_runs" on pipeline_runs
  for all using (auth.role() = 'service_role');

create policy "service_write_stages" on pipeline_stages
  for all using (auth.role() = 'service_role');

-- Enable Supabase Realtime for both tables
-- (Run this separately in Supabase dashboard if alter publication fails)
do $$
begin
  alter publication supabase_realtime add table pipeline_runs;
  alter publication supabase_realtime add table pipeline_stages;
exception
  when others then
    raise notice 'Could not add tables to realtime publication. Add them manually in Supabase dashboard.';
end;
$$;

-- Comments for documentation
comment on table pipeline_runs is 'Tracks each execution of a Marley bioinformatics pipeline';
comment on table pipeline_stages is 'Individual stage results within a pipeline run';
comment on column pipeline_runs.parameters is 'Frozen JSON of all input parameters for reproducibility';
comment on column pipeline_runs.tags is 'User-defined labels for filtering (e.g., dry-run, baseline)';
comment on column pipeline_stages.key_metrics is 'Stage-specific output metrics (e.g., proteins_found, duration)';
