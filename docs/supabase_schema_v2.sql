-- =============================================================================
-- Marley -- Supabase schema v2: vaccine construct tables (Module 06)
-- =============================================================================
-- Run this in the Supabase SQL Editor (or via psql) after the initial
-- candidates table has been created (see supabase_schema.sql).
-- =============================================================================

CREATE TABLE vaccine_constructs (
    construct_id        TEXT PRIMARY KEY,
    protein_sequence    TEXT NOT NULL,
    mrna_sequence       TEXT NOT NULL,
    signal_peptide      TEXT NOT NULL,
    adjuvant_name       TEXT NOT NULL,
    epitope_count       INTEGER NOT NULL,
    gc_content          REAL,
    molecular_weight    REAL,
    isoelectric_point   REAL,
    instability_index   REAL,
    gravy               REAL,
    vaxijen_score       REAL,
    allergenicity       TEXT,
    created_at          TIMESTAMPTZ DEFAULT now()
);

CREATE TABLE construct_epitopes (
    id                      SERIAL PRIMARY KEY,
    construct_id            TEXT REFERENCES vaccine_constructs(construct_id),
    peptide                 TEXT NOT NULL,
    epitope_type            TEXT NOT NULL,
    source_gene_id          TEXT,
    source_gene_name        TEXT,
    allele                  TEXT NOT NULL,
    ic50                    REAL NOT NULL,
    position_in_construct   INTEGER
);

-- Indexes for common queries
CREATE INDEX idx_constructs_created ON vaccine_constructs(created_at DESC);
CREATE INDEX idx_epitopes_construct ON construct_epitopes(construct_id);
