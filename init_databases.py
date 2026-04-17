"""
init_databases.py -- Initialize Brightway project and import databases

This script sets up a Brightway project and imports:
1. ecoinvent background database
2. Custom foreground databases (user-defined)

Usage:
    python init_databases.py

Before running, update the configuration section below with your own values.
"""
import os
import sys
import time
import bw2data as bd
import bw2io as bi

LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'init_db.log')


class Tee:
    def __init__(self, *files):
        self.files = files

    def write(self, data):
        for f in self.files:
            f.write(data)

    def flush(self):
        for f in self.files:
            f.flush()


log_f = open(LOG_FILE, 'w', encoding='utf-8')
sys.stdout = Tee(sys.stdout, log_f)
sys.stderr = Tee(sys.stderr, log_f)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ============================================================
# Configuration -- UPDATE THESE VALUES BEFORE RUNNING
# ============================================================
PROJECT_NAME = 'your_project_name'            # Brightway project name
ECOINVENT_VERSION = '3.11'                    # ecoinvent version (e.g. '3.11', '3.12')
ECOINVENT_SYSTEM_MODEL = 'consequential'       # System model: 'consequential' or 'cutoff'
ECOINVENT_USERNAME = os.getenv('ECOINVENT_USERNAME', 'your_username')  # ecoinvent username
ECOINVENT_PASSWORD = os.getenv('ECOINVENT_PASSWORD', 'your_password')  # ecoinvent password

# Custom databases to import after ecoinvent.
# Each entry: (db_name_in_bw, relative_path_to_xlsx, [list_of_dbs_to_match])
# The script will import them in the order listed.
CUSTOM_DATABASES = [
    # Example:
    # ('your_foreground_db', 'data/your_foreground.xlsx', []),
    # ('your_auxiliary_db',  'data/your_auxiliary.xlsx',  ['your_foreground_db']),
]


# ============================================================
# Step 1: Set project
# ============================================================
log("=" * 60)
log(f"Step 1: Setting project to '{PROJECT_NAME}'")
bd.projects.set_current(PROJECT_NAME)
log(f"  Current project: {bd.projects.current}")
log(f"  Existing databases: {list(bd.databases)}")


# ============================================================
# Step 2: Import ecoinvent if needed
# ============================================================
def find_ecoinvent_db():
    """Find the actual name of the ecoinvent DB matching configured version/model."""
    for name in bd.databases:
        name_lower = name.lower()
        if ('ecoinvent' in name_lower
                and ECOINVENT_SYSTEM_MODEL.lower() in name_lower
                and 'biosphere' not in name_lower):
            return name
    return None


log("")
log("=" * 60)
log(f"Step 2: Importing ecoinvent {ECOINVENT_VERSION} {ECOINVENT_SYSTEM_MODEL}")

existing_ei = find_ecoinvent_db()
if existing_ei:
    log(f"  ecoinvent already exists as '{existing_ei}', skipping.")
else:
    try:
        log(f"  Attempting to import ecoinvent {ECOINVENT_VERSION} {ECOINVENT_SYSTEM_MODEL}...")
        bi.import_ecoinvent_release(
            version=ECOINVENT_VERSION,
            system_model=ECOINVENT_SYSTEM_MODEL,
            username=ECOINVENT_USERNAME,
            password=ECOINVENT_PASSWORD,
        )
        log(f"  ecoinvent {ECOINVENT_VERSION} imported successfully.")
    except Exception as e:
        log(f"  FATAL: ecoinvent import failed: {e}")
        sys.exit(1)

log(f"  Databases after ecoinvent: {list(bd.databases)}")


# ============================================================
# Step 3: Import custom databases
# ============================================================
def find_ecoinvent_biosphere():
    """Find the actual name of the ecoinvent biosphere DB."""
    for name in bd.databases:
        if 'biosphere' in name.lower():
            return name
    return None


base_dir = os.path.dirname(os.path.abspath(__file__))
ei_consec = find_ecoinvent_db()
ei_bio = find_ecoinvent_biosphere()

if not ei_consec:
    log("WARNING: No ecoinvent database found. Custom DB imports may fail.")

log(f"  Using ecoinvent: {ei_consec}")
log(f"  Using biosphere: {ei_bio}")

for idx, (db_name, rel_path, match_dbs) in enumerate(CUSTOM_DATABASES, 1):
    log("")
    log("=" * 60)
    log(f"Step 3{chr(96 + idx)}: Importing '{db_name}'")

    if db_name in bd.databases:
        log(f"  '{db_name}' already exists, skipping.")
        continue

    abs_path = os.path.join(base_dir, rel_path)
    if not os.path.isfile(abs_path):
        log(f"  ERROR: File not found: {abs_path}")
        sys.exit(1)

    try:
        log(f"  Importing from: {abs_path}")
        ei = bi.ExcelImporter(abs_path)
        ei.apply_strategies()
        if ei_bio:
            ei.match_database(ei_bio)
        if ei_consec:
            ei.match_database(ei_consec)
        for other_db in match_dbs:
            if other_db in bd.databases:
                ei.match_database(other_db)
        ei.statistics()
        ei.write_project_parameters()
        ei.write_database(activate_parameters=True)
        log(f"  '{db_name}' imported successfully.")
    except Exception as e:
        log(f"  ERROR importing '{db_name}': {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


# ============================================================
# Final summary
# ============================================================
log("")
log("=" * 60)
log("FINAL SUMMARY")
log(f"  Project: {bd.projects.current}")
log(f"  Databases: {list(bd.databases)}")
for db_name in bd.databases:
    db = bd.Database(db_name)
    log(f"    {db_name}: {len(db)} activities")
log("Database initialization complete!")
log_f.close()
