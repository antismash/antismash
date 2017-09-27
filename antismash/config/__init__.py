# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import antismash.config.loader
import antismash.config.args

def update_config(values=None):
    return antismash.config.args.Config(values)

def get_config():
    return antismash.config.args.Config()

def destroy_config():
    antismash.config.args.Config().__dict__.clear()
