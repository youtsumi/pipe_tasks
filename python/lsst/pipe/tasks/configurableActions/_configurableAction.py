# This file is part of pipe_tasks.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
from __future__ import annotations

__all__ = ["ConfigurableAction"]

from typing import Any

from lsst.pex.config.config import Config


class ConfigurableAction(Config):
    """A `ConfigurableAction` is an interface that extends a
    `lsst.pex.config.Config` class to include a `__call__` method.

    This interfaced is designed to create an action that can be used at
    runtime with state that is determined during the configuration stage. A
    single action thus may be assigned multiple times, each with different
    configurations.

    In addition to allowing state to be set at configuration time these
    actions allow the state of execution to be recorded alongside all the
    other program configuration settings, allowing easy reproduction of
    results at some future time.

    This class is intended to be an interface only, and the `__call__` method
    is purely virtual. Subclasses that represent concrete actions must
    provide an override.
    """
    def __call__(self, *args, **kwargs) -> Any:
        raise NotImplementedError("This method should be overloaded in subclasses")
