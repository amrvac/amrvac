"""Velocity-estimation extension points for data-driven boundaries."""

from __future__ import print_function


class VelocityEstimator(object):
    """Protocol-like base class for future velocity estimators."""

    def estimate(self, previous_frame, current_frame):
        raise NotImplementedError


class Dave4VMEstimator(VelocityEstimator):
    """Placeholder for a future external DAVE4VM adapter.

    DAVE4VM is deliberately not vendored and not imported in V1. This class
    documents the future extension point while keeping non-DAVE workflows clean.
    """

    def __init__(self, *args, **kwargs):
        raise NotImplementedError(
            "DAVE4VM velocity estimation is not implemented in V1. "
            "Use B-only sequence preparation or provide a future adapter."
        )
