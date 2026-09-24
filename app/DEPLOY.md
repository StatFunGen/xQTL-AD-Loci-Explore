# AD xQTL Explorer — Deployment Guide

## Files
- `app.R`   — Shiny app (light/academic theme, full featured)
- `data.csv` — Processed data from Alexandre's unified table (4,796 rows)

## Option A: Deploy to shinyapps.io (recommended, free)

### 1. Install required packages (run once in R)
```r
install.packages(c("shiny", "DT", "dplyr", "readr", "stringr", "bslib"))
install.packages("rsconnect")
```

### 2. Connect your shinyapps.io account (run once)
```r
library(rsconnect)

# Get your token from: https://www.shinyapps.io/admin/#/tokens
rsconnect::setAccountInfo(
  name   = "YOUR_SHINYAPPS_USERNAME",
  token  = "YOUR_TOKEN",
  secret = "YOUR_SECRET"
)
```

### 3. Deploy
```r
rsconnect::deployApp(
  appDir  = "/path/to/shiny_app/",   # folder containing app.R and data.csv
  appName = "AD-xQTL-Explorer",
  account = "YOUR_SHINYAPPS_USERNAME"
)
```

Your app will be live at:
`https://YOUR_USERNAME.shinyapps.io/AD-xQTL-Explorer/`

### shinyapps.io free tier limits
- 5 apps
- 25 active hours/month (resets monthly)
- Fine for consortium sharing — upgrade to Starter ($9/mo) for unlimited hours

---

## Option B: Run locally
```r
shiny::runApp("/path/to/shiny_app/")
```

---

## Option C: GitHub Pages (HTML version, no server needed)
Use `AD_xQTL_academic.html` instead of the Shiny app:
1. Rename to `index.html`
2. Push to your `<user>/xqtl-paper` repo root (or a `docs/` folder)
3. Go to Settings → Pages → Source: main branch → `/` (or `/docs`)
4. Live at: `https://<user>.github.io/xqtl-paper/`

No server, no cost, no expiry. Best for a permanent consortium URL.

---

## Updating the data
When you get the full `.gz` flatten file and rerun your notebook:
1. Re-extract `data.csv` using the same Python script
2. For Shiny: re-deploy with `rsconnect::deployApp()`
3. For HTML: re-run the injection script and push the new `index.html`
