-- Presets: saved pipeline parameter configurations
-- Resolution order: user presets > team presets > system presets
-- System presets ship with the platform; users/teams can fork and customize
--
-- Depends on: 002_user_profiles.sql, 004_teams.sql

-- Table: presets
create table if not exists presets (
  id uuid primary key default gen_random_uuid(),
  pipeline text not null,
  name text not null,
  scope text not null default 'user'
    check (scope in ('system', 'team', 'user')),
  team_id uuid references teams(id) on delete cascade,
  created_by uuid references auth.users(id) on delete set null,
  display_name text not null default '',
  description text default '',
  recommended boolean not null default false,
  parameters jsonb not null default '{}',
  forked_from uuid references presets(id) on delete set null,
  created_at timestamptz not null default now(),
  updated_at timestamptz not null default now(),
  -- Constraints: scope determines which FK should be set
  constraint presets_scope_team_check
    check (
      (scope = 'system' and team_id is null)
      or (scope = 'team' and team_id is not null)
      or (scope = 'user' and team_id is null)
    ),
  -- Unique name per scope+pipeline+owner combination
  constraint presets_unique_name unique (pipeline, scope, team_id, created_by, name)
);

-- Auto-update updated_at
create trigger presets_updated_at
  before update on presets
  for each row execute function update_updated_at();

-- Indexes
create index if not exists idx_presets_pipeline on presets(pipeline);
create index if not exists idx_presets_scope on presets(scope);
create index if not exists idx_presets_team_id on presets(team_id) where team_id is not null;
create index if not exists idx_presets_created_by on presets(created_by) where created_by is not null;
create index if not exists idx_presets_pipeline_scope on presets(pipeline, scope);
create index if not exists idx_presets_forked_from on presets(forked_from) where forked_from is not null;

-- Row Level Security
alter table presets enable row level security;

-- System presets: visible to everyone (including anon for the wizard UI)
create policy "anyone_read_system_presets" on presets
  for select
  using (scope = 'system');

-- Team presets: visible to team members
create policy "members_read_team_presets" on presets
  for select to authenticated
  using (
    scope = 'team'
    and team_id = any(user_team_ids())
  );

-- User presets: visible only to their creator
create policy "users_read_own_presets" on presets
  for select to authenticated
  using (
    scope = 'user'
    and created_by = auth.uid()
  );

-- Users can create their own presets
create policy "users_insert_own_presets" on presets
  for insert to authenticated
  with check (
    created_by = auth.uid()
    and (
      scope = 'user'
      or (scope = 'team' and team_id = any(user_team_ids()))
    )
  );

-- Users can update their own presets
create policy "users_update_own_presets" on presets
  for update to authenticated
  using (created_by = auth.uid())
  with check (created_by = auth.uid());

-- Users can delete their own presets
create policy "users_delete_own_presets" on presets
  for delete to authenticated
  using (created_by = auth.uid());

-- Team admins can manage team presets
create policy "admins_manage_team_presets" on presets
  for all to authenticated
  using (
    scope = 'team'
    and team_id in (
      select tm.team_id from team_memberships tm
      where tm.user_id = auth.uid()
        and tm.role in ('owner', 'admin')
        and tm.accepted_at is not null
    )
  );

-- service_role full access (seeding system presets from backend)
create policy "service_all_presets" on presets
  for all using (auth.role() = 'service_role');

-- Comments
comment on table presets is 'Saved pipeline parameter configurations with user > team > system resolution';
comment on column presets.scope is 'Visibility: system (built-in), team (shared), user (personal)';
comment on column presets.parameters is 'Frozen JSON of pipeline parameters (same schema as pipeline_runs.parameters)';
comment on column presets.forked_from is 'Parent preset this was forked from (for lineage tracking)';
comment on column presets.recommended is 'Highlighted in the wizard UI as a suggested starting point';
comment on column presets.pipeline is 'Pipeline identifier this preset applies to (e.g., reverse_vaccinology)';
