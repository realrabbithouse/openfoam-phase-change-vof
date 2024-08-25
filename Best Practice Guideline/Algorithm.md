#***************************一些备注********************************#

1. 所有的mixture properties都应采用limitedAlpha1来进行相体积分数平均
2. 相变源项同样应该采用limitedAlpha1进行计算

#******************************************************************#

# interThermalPhaseChangeFoam Simulation Algorithm Summary

#******************************************************************#
interThermalPhaseChangeFoam的createFields.H和twoPhaseThermalMixture
都定义了mixture密度rho, 这样的操作似乎优点多余, 为了便于区分, 
在这里命名twoPhaseThermalMixture的mixture密度为rho.thermal
方程中用到的rho似乎都是rho
#******************************************************************#

WHILE runTime.run()
	1. 设置下一步时间步长
	2. runTime++
	3. 更新twoPhaseProperties: 通过limitedAlpha1更新mixture的nu_, lambda_, cp_, rho_ (rho.thermal)
	4. 更新muEffKistler [mixture动力粘度 + 湍流动力粘度]
	5. 更新phaseChangeModel: 更新Q_pc_
	(PS: 更新Q_pc_主要需要的是温度场T和相体积分数场alpha1, 显然这里采用了上一个时间步的T和alpha1)
	6. 求解相体积分数方程
		FOR iter < nAlphaSubCycles
			(PS: 求解compression velocity flux时代入的是上一个时间步的alpha1 field)
			FOR iter < nAlphaCorr
				通过MULES::explicitSolve求解alpha1 field
			END
		END
		更新interface properties (curvature K_ & interface normal nHatf_)
		更新mixture密度rho
	7. 根据interface properties和alpha1 field更新表面张力
	PIMPLE LOOP [nOuterCorrectors]
		1. 求解动量方程
		2. 求解压力方程
			WHILE pimple.correct() [nCorrectors]
				构造通量phi, 包含predicted velcoity flux, surface tension flux and gravitational force flux
				WHILE pimple.correctNonOrthogonal() [nNonOrthogonalCorrectors]
					求解p_rgh
					校正phi (仅在finalNonOrthogonalIter时校正, 此时phi是velocity volumetric flux)
				END
				校正速度U
				计算压力p, 然后计算压力p_rgh
			END
	8. turbulence->correct(): 通常是计算湍流的方程, i.e. k_ (kEqn), epsilon_ (epsEqn)
	9. 求解能量方程
		1. 计算等效导热系数和热扩散系数 (此时, mixture的导热系数和热扩散系数是通过alpha1更新过的)
		2. T.correctBoundaryConditions()
		3. 通过T计算H (质量平均而不是相体积分数平均)
		FOR iter < nEnergyLoops
			求解能量方程fvScalarMatrix EEqn
			通过H计算T
		END
	END
END

# 相变源项的显隐性

MULES::expilictSolve
	implicit = divU - mDot(1/rho1 - 1/rho2)
	explicit = mDot/rho1

Pressure Equation
	explicit: mDot(1/rho1 - 1/rho2)

Energy Equation
	explicit: Q_pc() = Lv*mDot



# phaseChangeHeatFoam Simulation Algorithm Summary

# 相变源项的显隐性

MULES::explicitSolve
	implicit = vDotAlphal[1] - vDotAlphal[0]
	explicit = alpha1*divU + vDotAlphal[0]
	mDot = (1 - alpha1)*mDotAlphal[0] + alpha1*mDotAlphal[1]

Pressure Equation
	implicit = vDotP[0] - vDotP[1]
	explicit = -(vDotP[0] - vDotP[1])*(pSat - rho*gh) // pSat是定值
	mDot = (mDotP[0] - mDotP[1])*(p - pSat)

Energy Equation
	implicit = (vDotT[0] - vDotT[1])
	explicit = (vDotT[0] - vDotT[1])*TSat // TSat由Clausius-Clapeyron方程确定
	mDot = (mDotT[0] - mDotT[1])*(T - TSat)



# interCondensatingEvaporatingFoam Simulation Algorithm Summary

WHILE runTime.run()
	1. 设置下一步时间步长
	2. runTime++
	WHILE pimple.loop() [nOuterCorrectors]
		1. 校正phase change model, 实际上什么也不做
		2. 求解相体积分数方程
			// ...
			更新mixture密度rho
		3. 更新interfaceProperties [interface normal nHatf_ & interface curvature K_]
		4. 求解连续性方程 [for mixture]
		5. 求解动量方程 // 动量方程当中用到了updated interfaceProperties和updated rho
		6. 求解能量方程 // 能量方程当中用到了updated mixture properties like Cp, kappa etc.
					  // 同时能量方程中的相变源项是基于updated limitedAlpha1的
		WHILE pimple.correct() [nCorrectors]
			1. 构造通量phiHByA [predicted velcoity flux, surface tension flux and gravitational force flux]
			// PS: gravitational force flux采用连续性方程更新后的mixture密度
			WHILE pimple.correctNonOrthogonal() [nNonOrthogonalCorrectors]
				1. 求解压力方程p_rghEqn
				2. 更新速度U和volumetric flux phi [仅在最后一个correctNonOrthogonal step]
			END
			2. 更新压力p, 并由p重新计算p_rgh
		END
		// PS: 压力方程中的源项采用了更新后的limitedAlpha1和上一次时间的温度T
		7. 求解湍流
	END
END



#*************************************isoAdvection*****************************************#

添加源项后需要修改的函数:
void advect();
void limitFluxes();
void boundFromAbove
(
    const scalarField& alpha1,
    surfaceScalarField& dVf, // 输入unbounded的timeIntegratedFlux
    DynamicList<label>& correctedFaces
);

PS:
(1) 关于VOF_和isoValue_的关系: 由VOF_计算一个cell的isoValue_, 并保证通过isoValue_来cut得到
的iso-surface恰好使得这个cell满足: subCellVolume_/mesh_.V()[celli] = VOF_
(2) fully submerged是指其状态为alpha1 = 1
(3) phi_ and dVf_ have same sign and dVf_ is the portion of phi_*dt that is liquid
(4) 注意cellIsBounded_与checkBounding_的区别, cellIsBounded_是checkBounding_的子集
 
##1 isoCutFace

Private Member Function:
// 计算一个多边形subFace的形心subFaceCentre_和面积subFaceArea_, subFaceArea_是vector类型
// 实际上subFace就是代表face中submerged的一部分; 计算的前提条件是subFacePoints_已知
void calcSubFaceCentreAndArea();

// 根据isoValue和face vertex VOF value计算cutEdge的系数, 并计数
// nFullySubmergedPoints_ (定义为f > isoValue),
// lastEdgeCut_定义在当f1 > isoValue > f2时, firstEdgeCut_定义在
// 当f1 < isoValue < f2时, 并且此时把f2定义为firstFullySubmergedPoint_
//  Returns the face status, where:
//  -1: face is fully below the isosurface
//   0: face is cut, i.e. has values larger and smaller than isoValue
//  +1: face is fully above the isosurface
label calcSubFace
(
    const scalar isoValue,
    const pointField& points, // mesh point
    const scalarField& f, // mesh point VOF value, f.size() = mesh.nPoints()
    const labelList& pLabels // face vertex index
)

// 把cutEdge对应的surfacePoints_和face所有的submerged point存入subFacePoints_中,
// 相当于保存了face的submerged部分对应的多边形
void subFacePoints
(
    const pointField& points, // mesh point
    const labelList& pLabels // face vertex index
)

// 通过cutEdge的系数lastEdgeCut_和firstEdgeCut_计算face的cutEdge得到的point,
// 计算通过firstFullySubmergedPoint_和nFullySubmergedPoints_来寻址需要的vertex,
// 结果保存在surfacePoints_中
void surfacePoints
(
    const pointField& points, // mesh point
    const labelList& pLabels // face vertex index
)

Public Member Function:
// 计算faceI的subFace相关信息
// i.e. subFacePoints_, subFaceArea_ and subFaceCenter_ etc.
label calcSubFace
(
    const label faceI,
    const scalar isoValue
)

// 这里的points和vertex VOF value f应该是特定的一个face的, 并不是通过全局的mesh vertex number来寻址
// 而是通过局部的face vertex number来寻址, i.e. 1, 2, 3..., 函数结果同上
label calcSubFace
(
    const pointField& points, // pointField of a face
    const scalarField& f, // vertex VOF value of a face
    const scalar isoValue
)

// 计算在time = f0时face-interface intersection line会运动到哪个位置, 并保存
// FIIL与face edge之间的交点到cutPoints中
void cutPoints
(
    const pointField& pts, // pointField of a concerned face
    const scalarField& f, // point arrival time, pTimes
    const scalar f0, // time
    DynamicList<point>& cutPoints
    // PS: 实际上用一个label构DynamicList时其size都初始化为0, 只是capacity会设置为label的值
)

// 已知初始时间和末时间的face-interface intersection lines FIIL和newFIIL, 计算FIIL和newFIIL
// 扫过的面积, 结果保存在系数alpha和beta中
void quadAreaCoeffs
(
    const DynamicList<point>& pf0, // FIIL
    const DynamicList<point>& pf1, // newFIIL
    scalar& alpha,
    scalar& beta
)

// timeIntegratedArea = integrate(t, t+dt, A(tau))
// where A(tau) is submerged area at time tau
Foam::scalar Foam::isoCutFace::timeIntegratedArea
(
    const pointField& fPts, // pointField of a face or a decomposed triangle
    const scalarField& pTimes,
    const scalar dt,
    const scalar magSf, // face area magnitude
    const scalar Un0 // velocity
)

// calculate volumetric face transport during dt given the isoFace
// data provided as input for facei
scalar timeIntegratedFaceFlux
(
	const label facei,
	const vector& x0, // interface center
	const vector& n0, // interface normal
	const scalar Un0, // velocity
	const scalar f0, // isoValue
	const scalar dt,
	const scalar phi,
	const scalar magSf
);

主要的疑惑:
(1) quadAreaCoeffs()中面积计算的原理以及为什么要那样积分不太明白
(2) timeIntegratedFaceFlux()中为什么要根据nShifts采用两种不同的计算策略不太明白

#********************************************************************************************#

##2 isoCutCell

Private Member Function:
// 已知isoCutFaceCentres_, isoCutFaceAreas_; fullySubFaces_; 以及isoFaceCentre_和
// isoFaceArea_, 计算subCellCentre_, subCellVolume_和VOF_
// isoCutFace: face is cut by isoface
// fullySubFace: face is fully submerged
// isoFace: interface between two phases
void calcSubCellCentreAndVolume()

// 已知isoFaceEdges_, 计算isoFaceCenter_和isoFaceArea_, 注意isoFaceArea_
// 指向永远是out of subCell
void calcIsoFaceCentreAndArea()

// 已知isoFaceEdges_, isoFaceArea_和isoFaceCenter_, 计算isoFacePoints_, 并且
// 按照角度来排序
void calcIsoFacePointsFromEdges()

Public Member Function:
// 已知isoValue, 计算celli所有face的isoCutFaces_相关的信息, 并统计fullySubFaces_;
// 如果cell is cut, 调用calcIsoFaceCentreAndArea(); 最后返回cellStatus_
label calcSubCell
(
    const label celli,
    const scalar isoValue
)

// 已知isofunction values at mesh points f_, alpha1 of celli, 计算celli应该的
// isoValue, 然后调用calcSubCell(celli, isoValue), 并返回cellStatus_
label vofCutCell
(
    const label celli,
    const scalar alpha1,
    const scalar tol,
    const label maxIter
)

#********************************************************************************************#

##3 isoAdvection

Private Member Function:

主要的疑惑:
(1) 为什么gradAlphaBasedNormal_时需要normalise?











