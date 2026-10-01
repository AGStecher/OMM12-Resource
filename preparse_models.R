library(xml2); library(dplyr)


# --- paste your full parse_gem_model() function here (with the get_attr fix) ---

# --- Parse a gapseq SBML genome-scale model into reactions + subsystems ---
parse_gem_model <- function(xml_path) {
  doc <- tryCatch(read_xml(xml_path), error = function(e) {
    message("GEM XML read error: ", e$message); NULL })
  if (is.null(doc)) return(NULL)
  
  ns   <- xml_ns(doc)
  core <- names(ns)[ns == "http://www.sbml.org/sbml/level3/version2/core"][1]
  fbc  <- names(ns)[ns == "http://www.sbml.org/sbml/level3/version1/fbc/version2"][1]
  grp  <- names(ns)[ns == "http://www.sbml.org/sbml/level3/version1/groups/version1"][1]
  if (is.na(core)) { message("Not an SBML core document"); return(NULL) }
  
  P <- function(tag, prefix = core) paste0(".//", prefix, ":", tag)
  
  # species id -> readable name (strip -c0/-e0/-p0 compartment suffix)
  sp <- xml_find_all(doc, P("species"))
  species_map <- setNames(sub("-[cep]0$", "", xml_attr(sp, "name")),
                          xml_attr(sp, "id"))
  
  # gene product id -> label
  gp <- xml_find_all(doc, P("geneProduct", fbc))
  gene_map <- setNames(xml_attr(gp, "label", ns = ns),
                       xml_attr(gp, "id",    ns = ns))
  
  # reactions
  rxn_nodes <- xml_find_all(doc, P("reaction"))
  fmt_side <- function(refs) {
    if (length(refs) == 0) return("")
    st   <- xml_attr(refs, "stoichiometry")
    sids <- xml_attr(refs, "species")
    nm   <- ifelse(sids %in% names(species_map), species_map[sids], sids)
    pref <- ifelse(is.na(st) | st %in% c("1", "1.0"), "", paste0(st, " "))
    paste(paste0(pref, nm), collapse = " + ")
  }
  rxn_list <- lapply(rxn_nodes, function(r) {
    tryCatch({
      rev <- isTRUE(xml_attr(r, "reversible") == "true")
      reactants <- fmt_side(xml_find_all(
        r, paste0("./", core, ":listOfReactants/", core, ":speciesReference")))
      products  <- fmt_side(xml_find_all(
        r, paste0("./", core, ":listOfProducts/",  core, ":speciesReference")))
      grefs <- xml_find_all(r, P("geneProductRef", fbc))
      gids  <- xml_attr(grefs, "geneProduct", ns = ns)
      glabs <- ifelse(gids %in% names(gene_map), gene_map[gids], gids)
      one <- function(x) if (length(x) == 0) NA_character_ else as.character(x)[1]
      
      data.frame(
        rid        = one(xml_attr(r, "id")),
        name       = one(xml_attr(r, "name")),
        reversible = isTRUE(rev),
        equation   = paste0(one(reactants), if (isTRUE(rev)) " <=> " else " => ", one(products)),
        genes      = if (length(glabs) == 0) "" else paste(sort(unique(glabs)), collapse = "; "),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message("GEM reaction parse failed at id=", xml_attr(r, "id"),
              " | error: ", conditionMessage(e))
      NULL
    })
  })
  rxn_list <- rxn_list[!vapply(rxn_list, is.null, logical(1))]
  reactions <- bind_rows(rxn_list)
  get_attr <- function(node, name, prefix) {
    v <- xml_attr(node, name, ns = ns)
    if (all(is.na(v))) v <- xml_attr(node, paste0(prefix, ":", name))  # explicit prefix
    if (all(is.na(v))) v <- xml_attr(node, name)                       # bare fallback
    v
  }
  
  # subsystems -> reaction members (long format)
  # subsystems -> members (long format) — read namespaced attributes explicitly
  sub_nodes <- if (!is.na(grp)) xml_find_all(doc, P("group", grp)) else list()
  sub_list <- lapply(sub_nodes, function(g) {
    sname <- get_attr(g, "name", grp)
    if (length(sname) == 0 || is.na(sname) || sname == "") {
      gid <- get_attr(g, "id", grp)
      sname <- sub("^subsys_", "", gid)   # fallback from groups:id
    }
    mems <- xml_find_all(g, P("member", grp))
    rids <- get_attr(mems, "idRef", grp)
    if (length(rids) == 0) return(NULL)
    data.frame(subsystem = sname, rid = rids, stringsAsFactors = FALSE)
  })
  subsystems <- bind_rows(sub_list)
  
  list(reactions = reactions, subsystems = subsystems)
}


# Parse every gapseq model found as data/{bacterium}_model.xml into the
# data/{bacterium}_model.rds the app reads. Names must match the app's
# bacteria_filenames exactly (including upper/lower case -- the server is Linux).
bacteria_filenames <- c(
  "acutalibacter_muris_kb18", "akkermansia_muciniphila_YL44",
  "bacteroides_caecimuris_I48", "bifidobacterium_animalis_YL2",
  "blautia_coccoides_YL58", "clostridium_innocuum_I46",
  "enterocloster_clostridioformis_YL32", "enterococcus_faecalis_KB1",
  "flavonifractor_plautii_YL31", "limosilactobacillus_reuteri_I49",
  "muribaculum_intestinales_YL27", "turicimonas_muris_YL45"
)

for (b in bacteria_filenames) {
  src <- file.path("data", paste0(b, "_model.xml"))
  out <- file.path("data", paste0(b, "_model.rds"))
  if (!file.exists(src)) next
  if (file.exists(out) && file.mtime(out) >= file.mtime(src)) {
    message("up to date: ", out); next
  }
  message("parsing ", src, " ...")
  t <- system.time(gem <- parse_gem_model(src))
  if (is.null(gem)) { message("  FAILED: ", src); next }
  saveRDS(gem, out)
  message(sprintf("  -> %s | %d rxns, %d subsystems | %.1fs",
                  out, nrow(gem$reactions),
                  length(unique(gem$subsystems$subsystem)), t["elapsed"]))
}
