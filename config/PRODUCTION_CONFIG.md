# Production Config

- Active frozen config file: `config/parameters_production.yaml`
- Freeze manifest: `config/production_freeze_2026_03_04.json`

Default config resolution order is now:
1. `BLACKCOTTON_CONFIG_PATH` (explicit override)
2. `config/parameters_production.yaml`
3. `config/parameters.yaml` (fallback)

This freeze was selected from the high-confidence anti-gibberish audit on 2026-03-04.
