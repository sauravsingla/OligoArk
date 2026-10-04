"""Train the transparent empirical policy baseline from reproducible observations."""

from oligoark.learning import EmpiricalPolicyModel, PolicyObservation
from oligoark.policy import ChannelProfile, CodecPolicy

lean = CodecPolicy(96, 8, 8, True, ("lean baseline",))
conservative = CodecPolicy(48, 24, 3, True, ("conservative baseline",))
observations = [
    PolicyObservation(ChannelProfile(), lean, True, 1200, 0.1),
    PolicyObservation(ChannelProfile(0.03, 0, 0, 0.10), lean, False, 1200, 0.1),
    PolicyObservation(ChannelProfile(0.03, 0, 0, 0.10), conservative, True, 2000, 0.2),
]
model = EmpiricalPolicyModel().fit(observations)
result = model.recommend(ChannelProfile(0.028, 0, 0, 0.09))
assert result.policy == conservative
print(result)
