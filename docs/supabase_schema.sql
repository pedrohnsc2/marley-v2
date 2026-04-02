-- =============================================================================
-- Marley -- Supabase schema for the candidates table
-- =============================================================================
-- Run this in the Supabase SQL Editor (or via psql) before starting the
-- pipeline for the first time.
-- =============================================================================

create table candidates (
  id                    uuid default gen_random_uuid() primary key,
  gene_id               text unique not null,
  gene_name             text,
  sequence              text,
  has_signal_peptide    boolean default false,
  conservation_score    float,
  immunogenicity_score  float,
  final_score           float,
  filters_passed        text[],
  status                text default 'pending',
  created_at            timestamp default now(),
  updated_at            timestamp default now()
);

-- Index for fast lookups by gene_id
create index idx_candidates_gene_id on candidates(gene_id);

-- Index for ranking queries
create index idx_candidates_final_score on candidates(final_score desc);
