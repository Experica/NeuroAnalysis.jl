using NeuroAnalysis.NACore
using Base.Test

v3 = Vec3()
av3 = [0.0,0.0,0.0]
@test v3 == convert(Vec3,av3)

