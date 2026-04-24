import math
from typing import Optional

from soundcalc.proxgaps.proxgaps_regime import ProximityGapsRegime


def _get_m_from_eta(sqrt_rate: float, eta: float) -> int:
    """
    Johnson-bound multiplicity m as a function of the gap η = 1 - √ρ - γ:

        m := max(⌈√ρ / (2η)⌉, 3)

    per BCHKS25 Theorem 1.5 / Theorem 4.2.

    The factor of 2 in the denominator is missing from the statement of BCHKS25 Theorem 4.2,
    which is a typo confirmed by the paper authors.
    """
    assert eta > 0, f"eta must be positive, got {eta}"
    return max(math.ceil(sqrt_rate / (2 * eta)), 3)


def _get_eta_from_m(sqrt_rate: float, m: int) -> float:
    """
    Inverse of `_get_m_from_eta`: η = √ρ / (2m). Drops the ⌈·⌉ and max-with-3 of
    the forward formula, so the caller who pinned m gets exactly the η that
    corresponds to their m.
    """
    assert m >= 1, f"m must be >= 1, got {m}"
    return sqrt_rate / (2 * m)


class JohnsonBoundRegime(ProximityGapsRegime):
    """
    Johnson Bound Regime (JBR).
    """
    def identifier(self) -> str:
        if self.explicit_m is not None:
            return f"JBR(m={self.explicit_m})"
        return "JBR"


    def __init__(
        self,
        field,
        gap_to_radius: Optional[float] = None,
        explicit_m: Optional[int] = None,
    ):
        super().__init__(field)

        assert not (gap_to_radius is not None and explicit_m is not None), (
            "gap_to_radius and explicit_m are mutually exclusive"
        )

        # Optional override for the Johnson-bound gap. If set, the proximity
        # parameter becomes: 1 - sqrt(rate) - gap_to_radius.
        # Primarily used by FRI-based systems that want fixed params.
        self.gap_to_radius = gap_to_radius
        # Optional override for the Johnson-bound multiplicity m. If set,
        # the proximity parameter and list size are derived directly from m
        # (rather than computing m from the gap). Used by protocols like
        # SWIRL that commit to a specific m.
        #
        # We allow m ≥ 1 (not just the m ≥ 3 from Thm 4.2): the lower bound
        # of 3 can be relaxed to 1, as noted in the optimization remark
        # right after Lemma 3.1 of BCHKS25.
        if explicit_m is not None:
            assert explicit_m >= 1, f"m must be >= 1, got {explicit_m}"
        self.explicit_m = explicit_m

    def get_proximity_parameter(self, rate: float, dimension: int) -> float:
        """
        The proximity parameter γ = 1 - √ρ - η from BCHKS25 Theorem 4.2.
        """
        sqrt_rate = math.sqrt(rate)

        # First let's pick our gap η
        if self.explicit_m:
            gap = _get_eta_from_m(sqrt_rate, self.explicit_m)
        elif self.gap_to_radius:
            gap = self.gap_to_radius
        else:
            # For large fields, use a tighter gap (closer to Johnson bound) for better
            # query-phase security. For smaller fields, use a more conservative gap.
            if self.field.F > 2**150:
                gap = sqrt_rate / 100
            else:
                gap = max(rate / 20, sqrt_rate / 100)

        assert(gap)
        return 1 - sqrt_rate - gap

    def get_max_list_size(self, rate: float, dimension: int) -> float:
        sqrt_rate = math.sqrt(rate)

        # With explicit m, the list size is (m + 0.5) / sqrt(rate).
        if self.explicit_m:
            return (self.explicit_m + 0.5) / sqrt_rate

        # Reed-Solomon codes are (1 - sqrt(rate) - gap, (2*gap*sqrt(rate))⁻¹)-list decodable.
        pp = self.get_proximity_parameter(rate, dimension)

        gap = 1 - sqrt_rate - pp
        assert gap > 0

        return 1.0 / (2 * gap * sqrt_rate)

    def get_m(self, rate: float, dimension: int) -> int:
        """
        Get m according to Theorem 4.2 of BCHKS25
        """
        if self.explicit_m:
            return self.explicit_m

        sqrt_rate = math.sqrt(rate)
        pp = self.get_proximity_parameter(rate, dimension)
        eta = 1 - sqrt_rate - pp
        return _get_m_from_eta(sqrt_rate, eta)

    def get_error_powers(self, rate: float, dimension: int, num_functions: int) -> float:
        return self.get_error_linear(rate, dimension) * (num_functions - 1)

    def get_error_linear(self, rate: float, dimension: int) -> float:
        """ Use Theorem 4.2 from BCHKS25 to compute the error"""

        sqrt_rate = math.sqrt(rate)

        pp = self.get_proximity_parameter(rate, dimension)
        m = self.get_m(rate, dimension)
        m_shifted = m + 0.5
        n = dimension / rate

        # Compute the first fraction
        numerator = (2 * m_shifted**5 + 3 * m_shifted * (pp * rate)) * n
        denominator = 3 * rate * sqrt_rate
        first_fraction = numerator / denominator

        # Now the second one
        second_fraction = m_shifted / sqrt_rate

        return (first_fraction + second_fraction) / self.field.F

    def get_error_multilinear(self, rate: float, dimension: int, num_functions: int) -> float:
        return self.get_error_linear(rate, dimension) * math.ceil(math.log2(num_functions))
