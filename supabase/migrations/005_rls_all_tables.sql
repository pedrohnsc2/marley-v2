-- RLS: team-based access across all tables
-- Replaces user-only policies from 003 with team-aware versions
--
-- Depends on: 004_teams.sql

-- Helper: returns array of team IDs for the current authenticated user
-- Used in RLS policies across multiple tables to avoid repeated subqueries
create or replace function user_team_ids()
returns uuid[] as $$
  select coalesce(
    array_agg(team_id),
    '{}'::uuid[]
  )
  from team_memberships
  where user_id = auth.uid() and accepted_at is not null
$$ language sql security definer stable;

-- ============================================================
-- pipeline_runs: upgrade to team-aware access
-- ============================================================

-- Drop the user-only policies from 003
drop policy if exists "authenticated_read_runs" on pipeline_runs;
drop policy if exists "authenticated_insert_runs" on pipeline_runs;

-- Read: own runs + team runs + legacy runs (no user_id)
create policy "team_read_runs" on pipeline_runs
  for select to authenticated
  using (
    user_id = auth.uid()
    or team_id = any(user_team_ids())
    or user_id is null  -- legacy runs before auth
  );

-- Insert: must be attributed to self, optionally within a team they belong to
create policy "team_insert_runs" on pipeline_runs
  for insert to authenticated
  with check (
    user_id = auth.uid()
    and (
      team_id is null
      or team_id = any(user_team_ids())
    )
  );

-- Update: own runs or runs within user's teams
create policy "team_update_runs" on pipeline_runs
  for update to authenticated
  using (
    user_id = auth.uid()
    or team_id = any(user_team_ids())
  )
  with check (
    user_id = auth.uid()
    or team_id = any(user_team_ids())
  );

-- service_role write policy from 001 remains unchanged

-- ============================================================
-- pipeline_stages: team-aware read through parent run
-- ============================================================

-- Drop the user-only policy from 003
drop policy if exists "authenticated_read_stages" on pipeline_stages;

-- Read: stages of runs the user can see (own + team + legacy)
create policy "team_read_stages" on pipeline_stages
  for select to authenticated
  using (
    exists (
      select 1 from pipeline_runs r
      where r.run_id = pipeline_stages.run_id
        and (
          r.user_id = auth.uid()
          or r.team_id = any(user_team_ids())
          or r.user_id is null
        )
    )
  );

-- service_role write policy from 001 remains unchanged

-- Comments
comment on function user_team_ids() is 'Returns team IDs for the current auth user (accepted memberships only). Used in RLS policies.';
