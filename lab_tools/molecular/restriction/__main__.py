"""Allow running as python -m lab_tools.molecular.restriction."""

import sys

from lab_tools.molecular.restriction.cli import main

sys.exit(main())
