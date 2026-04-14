-- Teams: multi-institution collaboration groups
-- Supports UFMG, Fiocruz, CTVacinas researcher collaboration
--
-- Depends on: 002_user_profiles.sql, 003_runs_add_user_id.sql

-- Table: teams
create table if not exists teams (
  id uuid primary key default gen_random_uuid(),
  name text not null,
  slug text not null unique,
  institution text default '',
  settings jsonb not null default '{}',
  created_at timestamptz not null default now(),
  updated_at timestamptz not null default now()
);

-- Auto-update updated_at
create trigger teams_updated_at
  before update on teams
  for each row execute function update_updated_at();

-- Indexes
create index if not exists idx_teams_slug on teams(slug);
create index if not exists idx_teams_institution on teams(institution) where institution is not null and institution != '';

-- Table: team_memberships
create table if not exists team_memberships (
  team_id uuid not null references teams(id) on delete cascade,
  user_id uuid not null references auth.users(id) on delete cascade,
  role text not null default 'researcher'
    check (role in ('owner', 'admin', 'researcher', 'viewer')),
  invited_by uuid references auth.users(id),
  accepted_at timestamptz,
  created_at timestamptz not null default now(),
  primary key (team_id, user_id)
);

-- Indexes for membership lookups
create index if not exists idx_team_memberships_user_id on team_memberships(user_id);
create index if not exists idx_team_memberships_team_id on team_memberships(team_id);

-- Now add the FK from pipeline_runs.team_id to teams
-- (column was added as plain uuid in 003, now we constrain it)
alter table pipeline_runs
  add constraint fk_pipeline_runs_team_id
  foreign key (team_id) references teams(id);

-- Row Level Security: teams
alter table teams enable row level security;

-- Members can see teams they belong to
create policy "members_read_teams" on teams
  for select to authenticated
  using (
    id in (
      select team_id from team_memberships
      where user_id = auth.uid() and accepted_at is not null
    )
  );

-- Owners and admins can update their team
create policy "admins_update_teams" on teams
  for update to authenticated
  using (
    id in (
      select team_id from team_memberships
      where user_id = auth.uid()
        and role in ('owner', 'admin')
        and accepted_at is not null
    )
  )
  with check (
    id in (
      select team_id from team_memberships
      where user_id = auth.uid()
        and role in ('owner', 'admin')
        and accepted_at is not null
    )
  );

-- Authenticated users can create teams (they become owner via app logic)
create policy "authenticated_insert_teams" on teams
  for insert to authenticated
  with check (true);

-- service_role full access
create policy "service_all_teams" on teams
  for all using (auth.role() = 'service_role');

-- Row Level Security: team_memberships
alter table team_memberships enable row level security;

-- Members can see memberships of their own teams
create policy "members_read_memberships" on team_memberships
  for select to authenticated
  using (
    team_id in (
      select tm.team_id from team_memberships tm
      where tm.user_id = auth.uid() and tm.accepted_at is not null
    )
  );

-- Users can see their own membership records (including pending invites)
create policy "users_read_own_memberships" on team_memberships
  for select to authenticated
  using (user_id = auth.uid());

-- Owners and admins can manage memberships
create policy "admins_manage_memberships" on team_memberships
  for all to authenticated
  using (
    team_id in (
      select tm.team_id from team_memberships tm
      where tm.user_id = auth.uid()
        and tm.role in ('owner', 'admin')
        and tm.accepted_at is not null
    )
  );

-- Users can update their own membership (accept invite)
create policy "users_accept_invite" on team_memberships
  for update to authenticated
  using (user_id = auth.uid())
  with check (user_id = auth.uid());

-- service_role full access
create policy "service_all_memberships" on team_memberships
  for all using (auth.role() = 'service_role');

-- Enable Realtime for teams (membership changes, settings updates)
do $$
begin
  alter publication supabase_realtime add table teams;
  alter publication supabase_realtime add table team_memberships;
exception
  when others then
    raise notice 'Could not add teams to realtime publication. Add them manually in Supabase dashboard.';
end;
$$;

-- Comments
comment on table teams is 'Collaboration groups spanning institutions (UFMG, Fiocruz, CTVacinas)';
comment on table team_memberships is 'Maps users to teams with role-based access';
comment on column teams.slug is 'URL-safe unique identifier (e.g., giunchetti-lab)';
comment on column teams.settings is 'Team-level configuration (default pipeline params, notification prefs)';
comment on column team_memberships.role is 'Team role: owner, admin, researcher, viewer';
comment on column team_memberships.accepted_at is 'Null until invite is accepted; used for access control';
