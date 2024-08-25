/*******************************************************************************\
                          When Rabbita Rules the World
                                A Rabbita Story
\*******************************************************************************/

带来方便: macro expansion || inline calculation and code


transformPoints -scale 0.001
splitMeshRegions -cellZones -overwrite

GFWList URL(Github): https://raw.githubusercontent.com/gfwlist/gfwlist/master/gfwlist.txt

Official mirror URLs:

Pagure: https://pagure.io/gfwlist/raw/master/f/gfwlist.txt

Repo.or.cz: http://repo.or.cz/gfwlist.git/blob_plain/HEAD:/gfwlist.txt

Bitbucket: https://bitbucket.org/gfwlist/gfwlist/raw/HEAD/gfwlist.txt

Gitlab: https://gitlab.com/gfwlist/gfwlist/raw/master/gfwlist.txt

TuxFamily: https://git.tuxfamily.org/gfwlist/gfwlist.git/plain/gfwlist.txt

Terminal使用代理：

export http_proxy=http://127.0.0.1:1088    有效

export ALL_PROXY=http://127.0.0.1:12333
export http_proxy=http://127.0.0.1:12333    有效
export https_proxy=http://127.0.0.1:12333

export ALL_PROXY=socks5://127.0.0.1:1080
export http_proxy=socks5://127.0.0.1:1080   http最好
export https_proxy=socks5://127.0.0.1:1080

g++编译C++代码：
g++ helloworld.cpp -o helloworld

// paraview label format
(i) %-#6.3e for exponential number
(ii) %-#6.3f for real number
(iii) %-#6.3g for integer number

// Meld is a visual diff and merge tool targeted at developers
// Meld helps you compare files, directories, and version controlled projects.
// It provides two- and three-way comparison of both files and directories, and
// has support for many popular version control systems
sudo apt-get install meld


if #include "dimensionedScalar.H" (always included)
then pMax_("pMax", dimPressure, great) can compile correctly


// B.C.
zeroGradientFvPatchScalarField::typeName
fixedValueFvPatchScalarField::typeName
extrapolatedCalculatedFvPatchScalarField::typeName

/**************************************************************************************************\
fvMesh: Mesh data needed to do the Finite Volume discretisation.

GeoMesh is a generic mesh wraper.
surfaceMesh: Mesh data needed to do the Finite Volume discretisation.
volMesh is derived from GeoMesh<fvMesh>
surfaceMesh is also derived from GeoMesh<fvMesh>

surfMesh: A surface mesh consisting of general polygon faces.
surfGeoMesh derived from GeoMesh<surfMesh>.

volScalarField: typedef GeometricField<scalar,fvPatchField,volMesh> volScalarField;
surfaceScalarField: typedef GeometricField<scalar, fvsPatchField, surfaceMesh> surfaceScalarField;

surfScalarField: typedef DimensionedField<scalar, surfGeoMesh> surfScalarField;
\**************************************************************************************************/

// alphaContactAngle
<patchName>
{
  // The essential entry "limit" controls the gradient of alpha1 on the wall:
  // none - Calculate the gradient from the contact-angle without limiter.
  // gradient - Limit the wall-gradient such that alpha1 remains bounded on the wall.
  // alpha - Bound the calculated alpha1 on the wall.
  // zeroGradient - Set the gradient of alpha1 to 0 on the wall, i.e. reproduce previous 
  //                behaviour, the pressure BCs can be left as before.
    type            constantAlphaContactAngle;
    theta0          90;
    limit           gradient;
    value           uniform 0;
}
<patchName>
{
    // advancing: increasing the volume, thetaA
    // receding:  decreasing the volume, thetaR
    // 一般: thetaA > thetaR
    type           dynamicAlphaContactAngle;
    uTheta         1;
    theta0         90;
    thetaA         70;
    thetaR         110;
    limit          gradient;
    value          uniform 0;
}
<patchName>
{
    type           dynamicAlphaContactAngle;
    uTheta         1;
    theta0         90;
    thetaAdv       110;
    thetaRec       70;
    limit          gradient;
    value          uniform 0;
}

heatingWall
{
  type            externalWallHeatFluxTemperature;
  mode            flux;
  q               uniform 1566586.397; // 100.41W
  // Ta              constant 300.0;
  // h               uniform 10.0;
  thicknessLayers (0.0009);
  kappaLayers     (401); // W/(m.K)
  kappaMethod     fluidThermo;
  value           $internalField;

  // type            fixedValue;
  // value           uniform 383;
}

coolingWall
{
    type            externalWallHeatFluxTemperature;
    mode            coefficient;
    Ta              constant 298.9;
    h               uniform 509.3;
    thicknessLayers (0.0009);
    kappaLayers     (401);
    kappaMethod     fluidThermo;
    value           $internalField;

    // type            fixedValue;
    // value           uniform 353;
}


/*******************************************************************************\
                                   常用命令
\*******************************************************************************/

nohup // 不挂断的运行命令
jobs // 查看当前终端的后台运行任务
ps // 查看瞬间进程的动态
kill // 结束进程
kill %jobnumber
kill %pidnumber
fg %jobnumber // 将后台任务调至前台

// to uncompress a file 
tar –xzvf filename.tar.gz
unzip filename.zip

// link
ln -s 源地址 目的地址

// man is short for manual, 一个命令的详细帮助信息
man cp

// install & uninstall
sudo dpkg -i // install
sudo dpkg -l // 查看安装的应用列表
sudo dpkg -r // uninstall

apt-get --purge remove <package> // 删除软件及其配置文件
apt-get autoremove <package> // 删除没用的依赖包
// 此时dpkg的列表中有“rc”状态的软件包，可以执行如下命令做最后清理：
dpkg -l |grep ^rc|awk '{print $2}' | sudo xargs dpkg -P

// find: d for directory; f for file
find $WM_PROJECT_DIR -type d -name “*fvPatch*”
find $WM_PROJECT_DIR -type f -name “*fvPatch*”

// find all the files that end with the Dict word in the tutorials directory
// -name for case sensitive search; -iname for non-case sensitive search
find $FOAM_TUTORIALS -name “*Dict”

// find files and scan them
find $FOAM_TUTORIALS -name controlDict | xargs grep 'probes' -sl

// scan a log file, '-c' means only display the number of matches
grep 'Max (alpha) = 1.' foamRun.log -c

// grep: r for recursive; n for line number
// find the string LES inside all the files within the directory $FOAM_SOLVERS
grep -r -n LES $FOAM_SOLVERS

// Find and scan files with file extension .C for the pattern
// lookupObject and count the hits
find $FOAM_SRC -name '*.C' | xargs grep 'lookupObject' | wc

// 批量替换satTempProperties -> saturationTemperatureProperties
grep -rl satTempProperties . | xargs sed -i "s/satTempProperties/saturationTemperatureProperties/g"
grep -rl constantLv . | xargs sed -i "s/constantLv/constant/g"
grep -rl constLatentHeat . | xargs sed -i "s/constLatentHeat/constant/g"
grep -rl empiricalParameter . | xargs sed -i "s/empiricalParameter/Lee/g"

// 将phaseChangeModels下文件中所有的Salpha替换为Salphal
sed -i "s/Salpha/Salphal/g" `grep Salpha -rl /phaseChangeModels`

// To find which tutorial files use the boundary condition "slip"
find $FOAM_TUTORIALS -type f | xargs grep -sl 'slip'

// diff is a command line tool that analyses two files 
// and prints a summary of the differences of those files.

// redirect
mpirun -np N icoFoam -parallel > foamRun.log
mpirun -np N icoFoam -parallel > /dev/ null

// To run cases with pyFoamPlotRunner.py, in the terminal type:
pyFoamPlotRunner.py icoFoam
// It is also possible to plot the information saved in the log file using PyFoam.
pyFoamPlotWatcher.py log.icoFoam

// chmod
chmod u=rwx,g=r,o=- file.txt
chmod u=rwx,g=r,o=r -R Application/ // 必须大写'R'

// gnuplot i.e.
set logscale y // set log scale in the y axis
gnuplot> plot ‘logs/p_0’ using 1:2 with lines
reset // reset the scales
// set the x range from 30 to 50 and plot tow files and set legend titles
gnuplot> plot [30:50][] ‘logs/Ux_0’ u 1:2 w l title ‘Ux’, ‘logs/Uy_0’ u 1:2 w l title ‘Uy’

// 死机时该怎么做？
1 直接 alt+f2 会弹出个输入框，输入 r 回车，这样会重启 gnome-shell 桌面环境

2 ctrl+f3 进入终端黑白屏环境
  top 一下，你会发现 gnome-shell cpu 100% 确认下眼神
  接着 kill pidNo
  然后系统还会自动重启桌面也就是 gnome-shell 稍等一会
  再按 ctrl+f1 切换会桌面环境登录即可

3 通用方法
  按住 alt+Prc Sc (SysRq) 然后按顺序按 r e i s u b 这样便会重启系统
  我的电脑: alt(一直按住)+fn+enter 然后按顺序按 r e i s u b

/******************************************************************************************\
                               Plot  information on the fly
             pyFoam: for details see p79 on OpenFOAM Tutorials-Wolf Dynamics
\******************************************************************************************/
 

/******************************************************************************************\
                                Tips & Tricks in OpenFOAM
\******************************************************************************************/       

// Dictionary Files General Features
// the lazy way
"(left|right|top)Wall"
{
	type fixedValue;
	value uniform (0 0 0);
}
".*Wall" 
{
	type fixedValue;
	value uniform (0 0 0);
}

// inline calculations using the directive #calc
X = 10.0; // declare variable
Y = 30.0;
Z  #calc  "$X*$Y - 12.0" // #eval

// multi-grading of a block, i.e.
blocks
(
    hex (0 1 2 3 4 5 6 7) (100 300 100)
    simpleGrading
    (
        1 // x-direction expansion ratio
        (
            (0.2 0.3 4) // 20% y-dir, 30% cells, expansion = 4
            (0.6 0.4 1) // 60% y-dir, 40% cells, expansion = 1
            (0.2 0.3 0.25) // 20% y-dir, 30% cells, expansion = 0.25 (1/4)
        )
        3 // z-direction expansion ratio
    )
);

// boundary faces ordering (expired)
The order in which the vertices are given must be such that,
looking from inside the block and starting with any vertex, the face
must be traversed in a clockwise direction to define the other vertices.

// Running in Parallel
decomposePar
mpirun -np 4 interFoam -parallel | tee log.interFoam
reconstructPar

// renumberMesh -Renumbers the cell list to reduce the bandwidth.
// autoRefineMesh -Refines cells near to a surface.
// refineHexMesh -Refines a hexahederal mesh by 2x2x2 cell splitting.

// Clear time directories
foamListTimes -rm

// 在编写OpenFOAM库时，通常将网格引用传递给基本模型的构造函数（实际上，如果没有引用网格，您将无能为力）

// Strictly speaking, the forAll macro just “works” on any container
// that has a size() member method.
forAll(list, i) // for (label i=0; i<(list).size(); ++i)

wmake - compiler seraching order for included files
1. the $WM_PROJECT_DIR/src/OpenFOAM/lnInclude directory;
2. a local lnInclude directory, i.e. newApp/lnInclude;
3. the local directory, i.e. newApp;
4. platform dependent paths set in files in the $WM_PROJECT_DIR/wmake/rules/-
   $WM_ARCH/ directory, e.g./usr/X11/include and $(MPICH_ARCH_PATH)/include;
5. other directories specified explicitly in the Make/options file with the -I option.


// Illustrate multi-region mesh generation process
ideasUnvToFoam fluid_copper_mesh.unv
transformPoints -scale 0.001
splitMeshRegions -cellZones -overwrite



/******************************************************************************************\
                                    Turbulence Model
\******************************************************************************************/

// Upper layer
transportModelIncompressibleTurbulenceModel == IncompressibleTurbulenceModel<transportModel>;
fluidThermoCompressibleTurbulenceModel == ThermoDiffusivity<CompressibleTurbulenceModel<fluidThermo>>;

// Lower layer
TurbulenceModel
<
    geometricOneField,
    geometricOneField,
    incompressibleTurbulenceModel,
    transportModel
>;
// Tables:
// laminarModel<IncompressibleTurbulenceModel<transportModel>>
// RASModel<IncompressibleTurbulenceModel<transportModel>>
// LESModel<IncompressibleTurbulenceModel<transportModel>>
TurbulenceModel
<
    geometricOneField,
    volScalarField,
    compressibleTurbulenceModel,
    fluidThermo
>；
// Tables:
// laminarModel<ThermalDiffusivity<CompressibleTurbulenceModel<fluidThermo>>>
// RASModel<EddyDiffusivity<ThermalDiffusivity<CompressibleTurbulenceModel<fluidThermo>>>>
// LESModel<EddyDiffusivity<ThermalDiffusivity<CompressibleTurbulenceModel<fluidThermo>>>>

// turbulenceModel.H
//- Return the turbulence kinetic energy
virtual tmp<volScalarField> k() const = 0;

//- Return the turbulence kinetic energy dissipation rate
virtual tmp<volScalarField> epsilon() const = 0;

//- Return the Reynolds stress tensor
virtual tmp<volSymmTensorField> R() const = 0;

// incompressibleTurbulenceModel.H
//- Return the effective stress tensor including the laminar stress
virtual tmp<volSymmTensorField> devReff() const = 0;

//- Return the source term for the momentum equation
virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const = 0;

// compressibleTurbulenceModel.H
//- Return the effective stress tensor including the laminar stress
virtual tmp<volSymmTensorField> devRhoReff() const = 0;

//- Return the source term for the momentum equation
virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const = 0;

// MULESCorr1
// 根据上一个时间步的dgdt和alpha1计算SuSp
// 组建alpha1Eqn, 更新alpha1并计算低阶迎风通量alphaPhi10
// 更新mixture (mixtrue.correct())

// nAlphaCorr
// 更新SuSp
// 根据alpha1显式计算没有限制器限制高阶通量talphaPhi1Un

// MULESCorr2
//     更新alpha1
//     计算corrected anti-diffusive flux
//     计算corrected higher-order flux
// MULES::explicitSolve
//     显式计算alpha1
//     显式计算corrected higher-order flux

// 更新mixture (mixtrue.correct())
// 计算rhoPhi

// PS:
// phir会随alpha1更新而更新
// volumetric flux phi在求解pEqn过后会更新

twoPhaseThermoMixture的接口
(1)
volScalarField& alpha1();
const volScalarField& alpha1() const;
const volScalarField& alpha1(mixture_.alpha1()); // const mixture_

(2)
virtual volScalarField& rho();
virtual tmp<volScalarField> rho() const;
const volScalarField& rho1 = mixture_.thermo1().rho(); // const mixture_

(3)
virtual const volScalarField& psi() const;
const volScalarField& psi1 = mixture.thermo1().psi(); // non-const mixture

(4)
virtual volScalarField& p();
virtual const volScalarField& p() const;
virtual volScalarField& T();
virtual const volScalarField& T() const;
volScalarField& p = mixture.p(); // non-const mixture
volScalarField& T = mixture.T(); // non-const mixture




