## Tips for All

1. Make sure that your p_rghFinal tolerance is well below the surfCellTol
Else you may be isoCutting non-surface cells which quadAreaCoeffs does not like


## VoFLibrary

1. 测试发现interFlow支持动网格, 参见$FOAM_RUN/VoFlibrary.verify/sloshingTank2D (mesh motion, no morphing change)
但不同的界面重构算法得到的结果明显不太一样, 结果显式isoAlpha的结果与iterIsoFoam的结果很相似
plicRDF的结果不太正确, 没有测试isoRDF, 建议最好动网格时用isoAlpha界面重构算法

2. damBreakWithObstacle (adaptive mesh refinement)测试同样发现, isoAlpha的结果与interIsoFoam的结果相似
但plicRDF同样算不下去, 日志显式出现warning, 且时间步极小


## isoAdvector

1. solve-continuity-problem branch is slightly different in function 'isoAdvection::advect()'


假设采用PISO循环
1. alphaEqn的作用
	更新相体积分数alpha1 & alpha2, 基于相体积分数更新mixture密度rho, 质量通量rhoPhi
	基于相体积分数更新mixture的其他物性参数, i.e. thermal diffusivity, mu, interface normal nHatf_, etc.

2. UEqn的作用
	动量预测, 然后更新比动能K

3. TEqn的作用
	根据mixture密度rho, 质量通量rhoPhi, 速度通量phi, 压力p, 比动能K, 连续性误差contErr等来更新mixture的温度T
	根据更新的温度T和压力p来校正mixture每一相的物性参数, i.e. psi, rho, mu, alpha, etc. (这说明受温度影响的物性与方程是完全耦合的)
	根据以上校正的每一相的物性参数来更新mixture的物性参数, 并更新interfaceProperties, 如果与温度有关的话

	PS: 如果TEqn在UEqn之后(OpenFOAM求解器中原本的处理方式), 那么动量预测后的比动能可能会带有很大的误差, 导致
	TEqn更新的温度也带有很大的误差!

4. pEqn做了什么?
	根据动量预测的速度等组建压力方程, 更新压力p_rgh & p, 然后更新速度U和速度通量phi
	最后, 根据压力的变化校正mixture每一相的密度, 然后更新mixture密度rho, 并更新比动能K

另一种方式(先完成速度压力的完全解耦后再求解TEqn):
1. main TEqn的位置调整
2. UEqn
	- K = 0.5magSqr(U);
3. pEqn
	+ rhoPhi = fvc::interpolate(rho)*phi
	+ contErr = (fvc::ddt(rho) + fvc::div(rhoPhi))()()

