-- API Keys: programmatic access for CI/CD and external tools
-- Stores hashed keys with scoped permissions per user or team
--
-- Depends on: 002_user_profiles.sql, 004_teams.sql

-- Table: api_keys
create table if not exists api_keys (
  id uuid primary key default gen_random_uuid(),
  user_id uuid not null references auth.users(id) on delete cascade,
  team_id uuid references teams(id) on delete cascade,
  key_hash text not null,
  key_prefix char(8) not null unique,
  name text not null default '',
  scopes text[] not null default '{run:create,run:read}',
  last_used_at timestamptz,
  expires_at timestamptz,
  revoked_at timestamptz,
  created_at timestamptz not null default now()
);

-- Indexes
create index if not exists idx_api_keys_user_id on api_keys(user_id);
create index if not exists idx_api_keys_team_id on api_keys(team_id) where team_id is not null;
create index if not exists idx_api_keys_key_prefix on api_keys(key_prefix);
create index if not exists idx_api_keys_key_hash on api_keys(key_hash);

-- Row Level Security
alter table api_keys enable row level security;

-- Users can read their own keys
create policy "users_read_own_keys" on api_keys
  for select to authenticated
  using (user_id = auth.uid());

-- Users can create keys for themselves
create policy "users_insert_own_keys" on api_keys
  for insert to authenticated
  with check (
    user_id = auth.uid()
    and (
      team_id is null
      or team_id = any(user_team_ids())
    )
  );

-- Users can update (revoke) their own keys
create policy "users_update_own_keys" on api_keys
  for update to authenticated
  using (user_id = auth.uid())
  with check (user_id = auth.uid());

-- Users can delete their own keys
create policy "users_delete_own_keys" on api_keys
  for delete to authenticated
  using (user_id = auth.uid());

-- Team admins can see keys scoped to their team
create policy "admins_read_team_keys" on api_keys
  for select to authenticated
  using (
    team_id in (
      select tm.team_id from team_memberships tm
      where tm.user_id = auth.uid()
        and tm.role in ('owner', 'admin')
        and tm.accepted_at is not null
    )
  );

-- service_role full access (key validation from backend)
create policy "service_all_keys" on api_keys
  for all using (auth.role() = 'service_role');

-- Comments
comment on table api_keys is 'Hashed API keys for programmatic pipeline access (CI/CD, scripts)';
comment on column api_keys.key_hash is 'SHA-256 hash of the full API key (plaintext never stored)';
comment on column api_keys.key_prefix is 'First 8 chars of the key for identification in UI (e.g., mrl_a1b2...)';
comment on column api_keys.scopes is 'Permitted operations: run:create, run:read, run:cancel, preset:write, etc.';
comment on column api_keys.revoked_at is 'Set when key is manually revoked; checked during validation';
comment on column api_keys.expires_at is 'Optional expiry; null means no expiry';
