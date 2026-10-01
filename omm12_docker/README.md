# OMM12 Docker / ShinyProxy deployment

Same structure as daisybio/namco's setup (https://github.com/daisybio/namco):
three containers on one network --

- **omm12** -- the R/Shiny app image (`app/Dockerfile`)
- **shinyproxy** -- starts a separate omm12 container for every visitor, so users
  run in parallel without sharing one R session (`shinyproxy/`)
- **nginx** -- reverse proxy in front of ShinyProxy, published on host port 3849
  (`nginx/nginx.conf`; namco uses 3848 on the same server)

The app is a pure data browser: it runs no bioinformatics tools, it only reads
the precomputed files in `data/`. All data is baked into the image (~290 MB of
app files after leaving out raw tool output).

## Layout

```
OMM12_website/
  .dockerignore          # decides what goes into the image: only app.R, data/, www/
  app.R  data/  www/
  omm12_docker/
    docker-compose.yml   # builds omm12 with the OMM12_website folder as build context
    app/Dockerfile
    shinyproxy/Dockerfile, application.yml
    nginx/nginx.conf
```

No copy/prepare step is needed any more: the image is built straight from the
project folder and `.dockerignore` leaves out raw tool output (`*_prokka/`,
`*_cog/`, `*_deeploc/`, `*.zip`), the E. coli Mt1B1 files, `www/data/`, and
everything outside app.R/data/www. (`prepare_build_context.sh` is obsolete.)

## Build and test locally (Windows: Docker Desktop, PowerShell)

```powershell
cd C:\Users\ge24biy\Documents\OMM12_website\omm12_docker
docker compose build          # first build takes ~10-20 min (R packages)
docker compose up -d
```

Then:

1. Open http://localhost:3849 -- you should land on the OMM12 app (first load
   takes a while, a new R container is started for you).
2. Open a second browser (or an incognito window) at the same time -- both should
   work independently. `docker ps` should show two `omm12:latest` containers.
3. Click through: OMM12 Resource grid, Bacteria Details -> Acutalibacter muris KB18
   (genome viewer incl. circular map, AMR section should show the vanR hit),
   Pan-genome, Relatedness Tree, Sequence Search.
4. Close the tabs; about a minute later the per-user containers should disappear
   from `docker ps`.
5. Memory per user: `docker stats` while a session is open -- useful number
   to send to the server admins.

Logs if something fails: `docker compose logs shinyproxy`, and the per-user
app logs in `omm12_docker/logs/`.

Stop everything: `docker compose down`.

Quick test of the app image alone (no ShinyProxy):
`docker run --rm -p 3838:3838 omm12:latest` -> http://localhost:3838

## Before handing over

- **R version**: `app/Dockerfile` has `ARG R_VERSION=4.4.1`. Set it to what
  `R.version.string` prints in your RStudio.
- **ShinyProxy admin password** in `shinyproxy/application.yml` (`CHANGE_ME`) --
  only protects ShinyProxy's admin page, the app itself is public.
- **/omm12 URL path**: this stack serves the app at `/` on port 3849. Ask the
  daisybio admins how `daisybio.ls.tum.de/omm12` will be routed to it (same as
  namco). If their front proxy forwards the `/omm12` prefix unchanged, add
  `server: servlet: context-path: /omm12` to `application.yml`.

## Changes vs. namco's config

- `container-wait-time: 90000` (namco: 10000) -- the app parses all 12 genomes
  and the ortholog tables at startup, which can take longer than 10 s.
- Network has a fixed name (`omm12_sp-net`), so `container-network` in
  application.yml matches regardless of the folder name it is deployed from.
- R packages install from rocker's dated package snapshot for the chosen R
  version (fixed versions, fast binary installs); the build fails early with a
  clear message if any package the app needs is missing.

## Troubleshooting

- `openjdk:11-jre` cannot be pulled -> in `shinyproxy/Dockerfile` use
  `FROM eclipse-temurin:11-jre` instead (the openjdk images are deprecated), and
  replace the `RUN wget ...` line with
  `ADD https://www.shinyproxy.io/downloads/shinyproxy-2.4.0.jar /opt/shinyproxy/shinyproxy.jar`
  (that image has no wget).
- ShinyProxy log says the Docker "client version ... is too old" (newer Docker
  engines dropped old API versions that ShinyProxy 2.4.0 uses) -> in
  `shinyproxy/Dockerfile` switch to a current ShinyProxy release, e.g.
  `FROM openanalytics/shinyproxy:3.1.1` (and remove the wget/jar lines).
