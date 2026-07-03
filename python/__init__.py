# Intentionally empty. The modules in this directory are standalone CLI scripts
# (run as `python <script>.py ...`) and are imported by the test suite via a
# sys.path shim (see tests/conftest.py), not as a package. Importing a submodule
# here would run import side effects during test collection.
