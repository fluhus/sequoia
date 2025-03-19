"""Local configuration."""

import json

_config = json.load(open('config.json'))

DATA_DIR = _config['dataDir']
"""General data directory."""

WS_DATA_DIR = _config['wsDataDir']
"""Workspace data directory."""

del _config
