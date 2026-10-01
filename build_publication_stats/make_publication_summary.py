#!/usr/bin/env python3
"""Aggregate the screened publication list into the small summary tables the
app's About page draws its charts from.

The per-publication list (data/omm12_publications.tsv) is built from a
Dimensions.ai free-version export plus Scopus citation counts. Their terms do
not allow redistributing that per-paper data, so it stays local (it is in
.gitignore and .dockerignore). Only these aggregated tables are published:

  data/omm12_pubstats_summary.tsv    key/value: totals, h-index, dates
  data/omm12_pubstats_years.tsv      publications + citations per year and group
  data/omm12_pubstats_countries.tsv  publications per author country (all / OMM12 users)
  data/omm12_pubstats_journals.tsv   top 10 peer-reviewed journals

Run after build_publication_stats.py (which calls this script automatically):
    python build_publication_stats/make_publication_summary.py
"""
import re
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"
GROUPS = ["Uses OMM12", "Other publications"]


def h_index(x):
    x = sorted([v for v in x if pd.notna(v)], reverse=True)
    return int(sum(v >= i + 1 for i, v in enumerate(x)))


def main():
    df = pd.read_csv(DATA / "omm12_publications.tsv", sep="\t", dtype=str, keep_default_na=False)
    df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
    df["citations_scopus"] = pd.to_numeric(df["citations_scopus"], errors="coerce")
    df["used"] = df["uses_omm12"].str.upper().eq("TRUE")
    df["group"] = df["used"].map({True: GROUPS[0], False: GROUPS[1]})

    info = {}
    p = DATA / "omm12_publications_info.tsv"
    if p.exists():
        i = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        info = dict(zip(i["key"], i["value"]))

    # publications by author country: each (publication, country) pair counts once
    rows = []
    for _, r in df.iterrows():
        iso = [s for s in re.split(r";\s*", r["countries_iso3"]) if s]
        nm = [s for s in re.split(r";\s*", r["countries"]) if s]
        for k, code in enumerate(iso):
            rows.append((code, nm[k] if k < len(nm) else code, r["used"]))
    links = pd.DataFrame(rows, columns=["iso3", "country", "used"])
    countries = (links.groupby(["iso3", "country"])
                 .agg(publications_all=("used", "size"), publications_uses_omm12=("used", "sum"))
                 .reset_index().sort_values(["publications_all", "country"], ascending=[False, True]))
    countries.to_csv(DATA / "omm12_pubstats_countries.tsv", sep="\t", index=False)

    # per year and group (every year x group present, zeros included)
    years = range(int(df["year"].min()), int(df["year"].max()) + 1)
    grid = pd.MultiIndex.from_product([years, GROUPS], names=["year", "group"])
    yearly = (df.groupby(["year", "group"])
              .agg(publications=("title", "size"), citations_scopus=("citations_scopus", "sum"))
              .reindex(grid, fill_value=0).reset_index())
    yearly["citations_scopus"] = yearly["citations_scopus"].astype(int)
    yearly.to_csv(DATA / "omm12_pubstats_years.tsv", sep="\t", index=False)

    # top 10 peer-reviewed journals (research articles and reviews)
    arts = df[(df["publication_type"] == "Article")
              & df["article_type"].isin(["Research article", "Review"]) & (df["journal"] != "")]
    totals = (arts.groupby("journal")
              .agg(journal_total=("title", "size"), journal_citations_scopus=("citations_scopus", "sum"))
              .reset_index().sort_values(["journal_total", "journal"], ascending=[False, True]).head(10))
    grid = pd.MultiIndex.from_product([totals["journal"], GROUPS], names=["journal", "group"])
    journals = (arts[arts["journal"].isin(totals["journal"])].groupby(["journal", "group"]).size()
                .reindex(grid, fill_value=0).rename("publications").reset_index()
                .merge(totals, on="journal"))
    journals["journal_citations_scopus"] = journals["journal_citations_scopus"].astype(int)
    journals.to_csv(DATA / "omm12_pubstats_journals.tsv", sep="\t", index=False)

    oa = df["open_access"].ne("") & df["open_access"].ne("Closed")
    summary = {
        "publications": len(df),
        "uses_omm12": int(df["used"].sum()),
        "year_min": int(df["year"].min()),
        "year_max": int(df["year"].max()),
        "countries": countries["iso3"].nunique(),
        "citations_scopus_total": int(df["citations_scopus"].sum()),
        "h_index_all": h_index(df["citations_scopus"]),
        "h_index_uses_omm12": h_index(df.loc[df["used"], "citations_scopus"]),
        "open_access_pct": round(100 * oa.mean()),
        "dimensions_export_date": info.get("dimensions_export_date", ""),
        "citations_date": info.get("citations_date", ""),
    }
    pd.DataFrame(list(summary.items()), columns=["key", "value"]).to_csv(
        DATA / "omm12_pubstats_summary.tsv", sep="\t", index=False)
    print("wrote data/omm12_pubstats_{summary,years,countries,journals}.tsv")


if __name__ == "__main__":
    main()
