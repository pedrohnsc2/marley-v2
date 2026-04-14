-- Runs: add user ownership and team association
-- Links pipeline_runs to authenticated users and (future) teams
--
-- Depends on: 001_pipeline_runs.sql, 002_user_profiles.sql

-- Add user_id: who launched this run
alter table pipeline_runs
  add column if not exists user_id uuid references auth.users(id);

-- Add team_id: which team context (populated once teams exist)
alter table pipeline_runs
  add column if not exists team_id uuid;

-- Indexes for ownership queries
create index if not exists idx_pipeline_runs_user_id on pipeline_runs(user_id);
create index if not exists idx_pipeline_runs_team_id on pipeline_runs(team_id) where team_id is not null;

-- Drop the old open-read policies from 001 so we can replace them
-- with user-scoped policies. Keep service_role write intact.
drop policy if exists "anon_read_runs" on pipeline_runs;
drop policy if exists "anon_read_stages" on pipeline_stages;

-- New RLS policies for pipeline_runs
-- Authenticated users see their own runs, team runs, or legacy runs (no user_id)
create policy "authenticated_read_runs" on pipeline_runs
  for select to authenticated
  using (
    user_id = auth.uid()
    or user_id is null  -- legacy runs before auth was added
  );

-- Authenticated users can create runs attributed to themselves
create policy "authenticated_insert_runs" on pipeline_runs
  for insert to authenticated
  with check (user_id = auth.uid());

-- service_role retains full access (Python worker backend)
-- (policy "service_write_runs" already exists from 001)

-- New RLS policies for pipeline_stages
-- Stages inherit visibility from their parent run
create policy "authenticated_read_stages" on pipeline_stages
  for select to authenticated
  using (
    exists (
      select 1 from pipeline_runs
      where pipeline_runs.run_id = pipeline_stages.run_id
        and (pipeline_runs.user_id = auth.uid() or pipeline_runs.user_id is null)
    )
  );

-- service_role retains full stage access from 001

-- Comments
comment on column pipeline_runs.user_id is 'Auth user who launched this run (null for legacy/pre-auth runs)';
comment on column pipeline_runs.team_id is 'Team context for this run (FK added in 004_teams migration)';
