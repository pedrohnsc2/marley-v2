-- User Profiles: researcher identity tied to Supabase Auth
-- Marley platform — multi-institution (UFMG, Fiocruz, CTVacinas)
--
-- Depends on: 001_pipeline_runs.sql (reuses update_updated_at trigger function)

-- Table: profiles
create table if not exists profiles (
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

-- Auto-update updated_at on row change (reuses function from 001)
create trigger profiles_updated_at
  before update on profiles
  for each row execute function update_updated_at();

-- Indexes
create index if not exists idx_profiles_email on profiles(email);
create index if not exists idx_profiles_orcid on profiles(orcid) where orcid is not null;
create index if not exists idx_profiles_institution on profiles(institution) where institution is not null and institution != '';

-- Auto-create profile on signup via Supabase Auth
create or replace function handle_new_user()
returns trigger as $$
begin
  insert into profiles (id, email, full_name)
  values (
    new.id,
    new.email,
    coalesce(new.raw_user_meta_data->>'full_name', '')
  );
  return new;
end;
$$ language plpgsql security definer;

create trigger on_auth_user_created
  after insert on auth.users
  for each row execute function handle_new_user();

-- Row Level Security
alter table profiles enable row level security;

-- Authenticated users can read any profile (needed for team member lists)
create policy "authenticated_read_profiles" on profiles
  for select to authenticated
  using (true);

-- Users can update only their own profile
create policy "users_update_own_profile" on profiles
  for update to authenticated
  using (id = auth.uid())
  with check (id = auth.uid());

-- Users can insert their own profile (fallback if trigger misses)
create policy "users_insert_own_profile" on profiles
  for insert to authenticated
  with check (id = auth.uid());

-- service_role has full access (backend operations)
create policy "service_all_profiles" on profiles
  for all using (auth.role() = 'service_role');

-- Comments
comment on table profiles is 'Researcher profiles linked to Supabase Auth users';
comment on column profiles.role is 'Platform role: admin, pi (principal investigator), researcher, viewer';
comment on column profiles.orcid is 'ORCID identifier for academic attribution';
comment on column profiles.institution is 'Primary institution (e.g., UFMG, Fiocruz, CTVacinas)';
