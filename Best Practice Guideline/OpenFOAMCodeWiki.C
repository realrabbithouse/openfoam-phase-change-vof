// fvPatchField & fvsPatchField
fvPatchField gives boundary conditions for volField<Type>
fvsPatchField describes boundary conditions for surfaceField<Type>

typedef GeometricField<scalar, fvPatchField, volMesh> volScalarField;
typedef GeometricField<scalar, fvsPatchField, surfaceMesh> surfaceScalarField

// adjusts the inlet and outlet fluxes to obey continuity, which is necessary for
// creating a well-posed problem where a solution for pressure exists.
// adjustPhi 输入的 face flux 'phi' is relative to mesh motion and MRF!
adjustPhi(phi, U, p);

// The ddtPhiCorr term accounts for the divergence of the face velocity field by taking out the 
// difference between the interpolated velocity and the flux.
fvc::ddtPhiCorr(rUA, U, phi);

extern const dimensionSet dimless;
extern const dimensionSet dimMass;
extern const dimensionSet dimLength;
extern const dimensionSet dimTime;
extern const dimensionSet dimTemperature;
extern const dimensionSet dimMoles;
extern const dimensionSet dimCurrent;
extern const dimensionSet dimLuminousIntensity;

extern const dimensionSet dimArea;
extern const dimensionSet dimVolume;
extern const dimensionSet dimVol;

extern const dimensionSet dimDensity;
extern const dimensionSet dimForce;
extern const dimensionSet dimEnergy;
extern const dimensionSet dimPower;

extern const dimensionSet dimVelocity;
extern const dimensionSet dimAcceleration;
extern const dimensionSet dimPressure;
extern const dimensionSet dimCompressibility;
extern const dimensionSet dimGasConstant;
extern const dimensionSet dimSpecificHeatCapacity;
extern const dimensionSet dimViscosity;
extern const dimensionSet dimDynamicViscosity;

/******************************************************************************************\
                                runTimeSelection mechanism
\******************************************************************************************/

declareRunTimeSelectionTable(autoPtr,baseType,argNames,argList,parList); // 1st
/* 
   argNames - this is a small name you gave to the parameter list (i.e. MrConstructor below).
   argList - this is the full parameter list with modifiers, types and names - enclose in ()'s.
   parList - this is the parameter list with names only. 
*/
// In Base.H
// 构造函数指针, 返回类型是Base
typedef Base (*MrConstructorPtr)( const dictionary& dict );
// or, i.e.
typedef Base (*MrsConstructorPtr)
        (
            const label& place,
            const scalar& magnitude,
            const Istream& dataFlow
        );

// Showing only MrConstructor from now on
// (word -> 构造函数指针)的hashTable
typedef HashTable<MrConstructorPtr, word, string::hash>
    MrConstructorTable;
static MrConstructorTable* MrConstructorTablePtr_; // hashTable的指针

// prototypes for the hash table creator and destroyer functions
static void constructMrConstructorTables();
static void destroyMrConstructorTables();

// the full definition of the AddToTable subclass
template<class DerivedType>
class addMrConstructorToTable
{
public:
 
    static Base New(const dictionary& dict)
    {
        return Base( new DerivedType(dict) );
    }
 
    addMrConstructorToTable
    (
        const word& lookup = DerivedType::typeName
    )
    {
        constructMrConstructorTables();
        // hashTable的insert(key, value)方法
        MrConstructorTablePtr_->insert(lookup, New);
    }
 
    ~addMrConstructorToTable()
    {
        destroyMrConstructorTables();
    }
};


defineRunTimeSelectionTable(baseType,argNames); // 2nd

// In Base.C
void Base::constructMrConstructorTables()
{
    // static variable, when initialized, it will exist untill program termination
    static bool constructed = false;
    
    if (!constructed)
    {
        // 空的hashTable
        Base::MrConstructorTablePtr_ = new Base::MrConstructorTable;
 
        constructed = true;
    }
}
 
void Base::destroyMrConstructorTables()
{
    if (Base::MrConstructorTablePtr_)
    {
        delete Base::MrConstructorTablePtr_;
        Base::MrConstructorTablePtr_ = NULL;
    }
}


addToRunTimeSelectionTable(baseType,thisType,argNames); // 3th

// In Derived.C
Base::addMrConstructorToTable<DerivedType>
    addDerivedTypeMrConstructorToBaseTable_;


// Selector - the virtual constructor function // 4th, no macro

// In Base.H
static tmp<Base> New(const dictionary&);

// In Base.C(or BaseNew.C)
Foam::tmp<Foam::Base> Foam::Base::New(const dictionary& dict)
{
    // omitting error catching and debug statements
 
    word DerivedType(dict.lookup("type"));
 
    typename MrConstructorTable::iterator cstrIter
        = MrConstructorTablePtr_->find(DerivedType);
 
    return cstrIter()(dict);
}


/******************************************************************************************\
                        	TypeName & defineTypeNameAndDebug
\******************************************************************************************/

// TypeName("transportModel") Macro

 #define TypeName(TypeNameString)                                               \
     ClassName(TypeNameString);                                                 \
     virtual const word& type() const { return typeName; }

 #define ClassName(TypeNameString)                                              \
     ClassNameNoDebug(TypeNameString);                                          \
     static int debug

 #define ClassNameNoDebug(TypeNameString)                                       \
     static const char* typeName_() { return TypeNameString; }                  \
     static const ::Foam::word typeName

// TypeName("transportModel") Macro unfold
static const char* typeName_() { return "transportModel"; }
static const ::Foam::word typeName;
static in debug;
virtual const word& type() const { return typeName; }



// defineTypeNameAndDebug(transportModel, 0) Macro
//- Define the typeName and debug information
#define defineTypeNameAndDebug(Type, DebugSwitch)                              \
    defineTypeName(Type);                                                      \
    defineDebugSwitch(Type, DebugSwitch)

//- Define the typeName
#define defineTypeName(Type)                                                   \
    defineTypeNameWithName(Type, Type::typeName_())

//- Define the typeName, with alternative lookup as \a Name
#define defineTypeNameWithName(Type, Name)                                     \
    const ::Foam::word Type::typeName(Name)

//- Define the debug information
#define defineDebugSwitch(Type, DebugSwitch)                                   \
    defineDebugSwitchWithName(Type, Type::typeName_(), DebugSwitch);           \
    registerDebugSwitchWithName(Type, Type, Type::typeName_())

//- Define the debug information, lookup as \a Name
#define defineDebugSwitchWithName(Type, Name, DebugSwitch)                     \
    int Type::debug(::Foam::debug::debugSwitch(Name, DebugSwitch))

//- Define the debug information, lookup as \a Name
#define registerDebugSwitchWithName(Type,Tag,Name)                             \
    class add##Tag##ToDebug                                                    \
    :                                                                          \
        public ::Foam::simpleRegIOobject                                       \
    {                                                                          \
    public:                                                                    \
        add##Tag##ToDebug(const char* name)                                    \
        :                                                                      \
            ::Foam::simpleRegIOobject(Foam::debug::addDebugObject, name)       \
        {}                                                                     \
        virtual ~add##Tag##ToDebug()                                           \
        {}                                                                     \
        virtual void readData(Foam::Istream& is)                               \
        {                                                                      \
            Type::debug = readLabel(is);                                       \
        }                                                                      \
        virtual void writeData(Foam::Ostream& os) const                        \
        {                                                                      \
            os << Type::debug;                                                 \
        }                                                                      \
    };                                                                         \
    add##Tag##ToDebug add##Tag##ToDebug_(Name)

// defineTypeNameAndDebug(transportModel, 0) Macro unfold
// transportModel::typyName_() = "transportModel"
const ::Foam::word transportModel::typeName("transportModel")
// defineDebugSwitchWithName(transportModel,"transportModel",0);
int transportModel::debug(::Foam::debug::debugSwitch("transportModel",0))
// registerDebugSwitchWithName(transportModel,transportModel,"transportModel");
class addtransportModelToDebug_
:
    public ::Foam::simpleRegIOobject
{
public:
    addtransportModelToDebug(const char* name)
    :
        ::Foam::simpleRegIOobject(Foam::debug::addDebugObject, name)
    {}
    virtual ~addtransportModelToDebug()
    {}
    virtual void readData(Foam::Istream& is)
    {
        transportModel::debug = readLabel(is);
    }
    virtual void writeData(Foam::Ostream& os) const
    {
        os << transportModel::debug;
    }
};
addtransportModelToDebug addtransportModelToDebug_("transportModel");

/******************************************************************************************\

\******************************************************************************************/

tmp<volScalarField> tVSF
(
    volScalarField::New
    (
        "VSF",
        mesh,
        dimensionedScalar(/.../),
        patchFieldType // 默认实参: PatchField<Type>::calculatedType()
    )
);  
// 等价操作
tmp<volScalarField> tVSF
(
    new volScalarField
    (
        IOobject
        (
            "VSF",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(/.../),
        patchFieldType
    )
)


// Correct Uf if the mesh is moving
// 这样做的实际意义是使得Uf在面法向的投影的大小等于phi/mesh.magSf()
fvc::correctUf(Uf, U, phi);
const fvMesh& mesh = U.mesh();

if (mesh.dynamic())
{
   Uf() = fvc::interpolate(U);
   surfaceVectorField n(mesh.Sf()/mesh.magSf());
   Uf() += n*(phi/mesh.magSf() - (n & Uf()));
}

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);
if (phi.mesh().moving())
{
    phi -= fvc::meshPhi(U);
}

// Return the given relative flux in absolute form
// 返回临时的absolute flux, 实际意义并不像fvc::makeRelative(phi, U)
fvc::absolute(phi, U);
if (tphi().mesh().moving())
{
   return tphi + fvc::meshPhi(U);
}
else
{
   return tmp<surfaceScalarField>(tphi, true);
}


mesh.dynamic(): true for all dynamicFvMesh except staticFvMesh
mesh.moving(): Is mesh moving?
mesh.topoChanging(): Is mesh topology changing?
mesh.changing(): Is mesh changing (topology changing and/or moving)?


// IO 之 直接写入
if (mesh_.time().outputTime())
{
    mDotC.write();
    mDotE.write();
}

相同的表达式fvMatrix<Type>和volField<Type>的单位是有区别的, i.e.
volScalarField q1 = coeff*T; // 单位 = unit
fvScalarMatrix q2 = fvm::(coeff, T); // 单位 = unit*m^3
这可能与有限体积法离散有关系

// isoAdvection access interface normal
// interface normal is updated after interface reconstruction step
isoAdvector.surf().reconstruct();
isoAdvector.surf().normal();

/******************************************************************************************\
                                   thermoPhysical Model
\******************************************************************************************/
// thermophysicalProperties
thermoType
{
    // thermalphysical model
    type            hePsiThermo;

    // The mixture specifies the mixture composition. The option typically used for thermophysical 
    // models without reactions is pureMixture, which represents a mixture with fixed composition.
    mixture         pureMixture;

    // transport model: evaluating dynamic viscosity, thermal conductivity and thermal diffusivity.
    transport       const;

    // thermodynamic model: evaluating the specific heat Cp from which other properties are derived.
    // maybe specifiy a heat of fusion coefficient Hf (i.e. hConst, eConst)
    thermo          hConst;

    // equation of state: evaluating density.
    equationOfState perfectGas;

    // There is currently only one option for the specie model which specifies the composition of
    // each constituent.
    // nMoles: this entry is only used for combustion modelling.
    // molWeight: grams per mole of specie.
    specie          specie;

    // we refer to absolute energy where heat of formation is included, i.e. absoluteEnthalpy.
    energy          sensibleEnthalpy; // sensibleInternalEnergy, absoluteEnthalpy
}
thermoTypeName == "heRhoThermo<pureMixture<const<hConst<perfectGas<specie>>,sensibleEnthalpy>>"

// #define makeThermo(BaseThermo,Cthermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie)
// i.e.
makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);
// 1) typedef
// #define typedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie)
typedef constTransport<speies::thermo<hConstThermo<perfectGas<specie>>,sensibleInternalEnergy>>
constTransportsensibleInternalEnergyhConstThermoperfectGasspecie;
// Transport##Type##Thermo##EqnOfState##Specie

// 2) typedef & defineTemplateTypeNameAndDebugWithName
// #define defineThermoPhysicsThermo(BaseThermo,CThermo,Mixture,ThermoPhys)
// ThermoPhys == constTransportsensibleInternalEnergyhConstThermoperfectGasspecie
typedef heRhoThermo
        <
            rhoThermo,
            pureMixture<constTransportsensibleInternalEnergyhConstThermoperfectGasspecie>
        > heRhoThermopureMixtureconstTransportsensibleInternalEnergyhConstThermoperfectGasspecie;
          // Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie
defineTemplateTypeNameAndDebugWithName
(
    heRhoThermopureMixtureconstTransportsensibleInternalEnergyhConstThermoperfectGasspecie,
    "heRhoThermo<pureMixture<const<hConst<perfectGas<specie>>,sensibleInternalEnergy>>>",
    0
) // Specie::typeName_() based on macro ClassName("specie"), which defines a function typeName_()
  // Specie::typeName_() -> "specie"

// 3) #define addThermoPhysicsThermo(BaseThermo,CThermoMixtureThermoPhys)
// CThermoMixtureThermoPhys == 
// heRhoThermopureMixtureconstTransportsensibleInternalEnergyhConstThermoperfectGasspecie
addToRunTimeSelectionTable
(
    rhoThermo,
    heRhoThermopureMixtureconstTransportsensibleInternalEnergyhConstThermoperfectGasspecie,
    fvMesh
); 
// addToRunTimeSelectionTable的key由该类的typeName确定，
// 类heRhoThermopureMixtureconstTransportsensibleInternalEnergyhConstThermoperfectGasspecie的key是
// heRhoThermopureMixtureconstTransportsensibleInternalEnergyhConstThermoperfectGasspecie::typeName
// 实际上是"heRhoThermo<pureMixture<const<hConst<perfectGas<specie>>,sensibleInternalEnergy>>"，在宏
// defineTemplateTypeNameAndDebugWithName当中定义

// 同时需要注意模板类的typeName和模板类实例化后的类的typeName的区别，如：
TypeName("heRhoThermo"); // 声明了模板类的typeName但不定义
defineTemplateTypeNameAndDebugWithName
(
    heRhoThermopureMixtureconstTransportsensibleInternalEnergyhConstThermoperfectGasspecie,
    "heRhoThermo<pureMixture<const<hConst<perfectGas<specie>>,sensibleInternalEnergy>>",
    0
) // 定义了模板类实例化后的类的typeName

// 总结：上述makeThermo的调用，向rhoThermo类中声明的hashTable中添加一组元素，其key为
// "heRhoThermo<pureMixture<const<hConst<perfectGas<specie>>,sensibleInternalEnergy>>>"
// value对应的函数返回的是类
// heRhoThermo                                                                 
// <                                                                        
//     rhoThermo,                                                          
//     pureMixture
//     <
//         constTransport<species::thermo<hConstThermo<perfectGas<specie>>,sensibleInternalEnergy>>
//     >                 
// >
// 的对象.

// Fundamental properties in mass specific (thermo.H) | species::thermo
// 这些函数实际上在species::thermo的原则类当中定义的

            Heat capacity at constant pressure [J/kg/K]
            inline scalar Cp(const scalar p, const scalar T) const;

            Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            Heat capacity at constant volume [J/kg/K]
            inline scalar Cv(const scalar p, const scalar T) const;

            Sensible internal energy [J/kg]
            inline scalar Es(const scalar p, const scalar T) const;

            Absolute internal energy [J/kg]
            inline scalar Ea(const scalar p, const scalar T) const;

            Entropy [J/kg/K]
            inline scalar S(const scalar p, const scalar T) const;

// Derived properties (mass specific)

            //- Heat capacity at constant pressure/volume [J/kg/K]
            inline scalar Cpv(const scalar p, const scalar T) const;

            //- Gamma = Cp/Cv []
            inline scalar gamma(const scalar p, const scalar T) const;

            //- Ratio of heat capacity at constant pressure to that at
            //  constant pressure/volume []
            inline scalar CpByCpv(const scalar p, const scalar T) const;

            //- Enthalpy/Internal energy [J/kg]
            inline scalar HE(const scalar p, const scalar T) const;

            //- Gibbs free energy [J/kg]
            inline scalar G(const scalar p, const scalar T) const;

            //- Helmholtz free energy [J/kg]
            inline scalar A(const scalar p, const scalar T) const;

// Energy->temperature inversion functions

            //- Temperature from enthalpy or internal energy
            //  given an initial temperature T0
            inline scalar THE
            (
                const scalar H,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from sensible enthalpy given an initial T0
            inline scalar THs // THa for temperature from absolute enthalpy
            (
                const scalar Hs,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
            inline scalar TEs // TEa for temperature from absolute internal energy
            (
                const scalar E,
                const scalar p,
                const scalar T0
            ) const;

// what's inside?
class specie
{
    word name_; // this name_ is a dictionary name, not a phaseName!!!
    //- Number of moles of this component in the mixture
    scalar Y_;
    scalar molWeight_; // [kg/kmol]
    // key method: specific gas constant R() & molecular weight W()
} // Cp - Cv = R/W
template class perfaceGas<class Specie> // equationOfState
{
    // key method: density rho(p, T) & compressibility psi(p, T) 
}
template class hConstThermo<class EquationOFState> // thermo
{
    scalar Cp_; // or scalar Cv_
    scalar Hf_; // heat of formation 
}
template<class Thermo, template<class> class Type> class thermo // species::thermo
{
    static const scalar tol_;
    static const int maxIter_;
    // key method: calculate T based on energy type, i.e.
    // for sensibleEnthalpy, it will calculate via THs(...) method
}
template<class Thermo> class constTransport
{
    scalar mu_;
    scalar rPr_; // 1/Pr
}
class basicThermo // 继承自IOdictionary
{
    const word& phaseName_; // 引用类型！！！
    volScalarField& p_;
    volScalarField T_; // T_的名字是"T.phaseName_"

    volScalarField alpha_; // liminar thermal diffusivity
    // PS: alpha_的单位是[kg/m/s], 注意与specific thermal diffusivity的区别[m^2/s]

    Switch dpdt_; // for energy equation based on enthalpy
}
class rhoThermo // 继承自fluidThermo
{
    volScalarField rho_;
    volScalarField psi_;
    volScalarField mu_;
}
class psiThermo
{
    volScalarField psi_;
    volScalarField mu_;
}
template<class BasicThermo, class MixtureType> class heThermo
{
    volScalarField he_; // he_的名字是i.e. "h.phaseName_"
}

PS:
(1) pureMixture什么也不做，只是作为一个ThermoType的wraper，i.e. wraper for specialized template class
    constTransport<species::thermo<hConstThermo<perfectGas<specie>>,sensibleInternalEnergy>>
(2) pureMixture说明ThermoType的properties从subDict "mixture"字典中读取的
(3) 完全类(i.e. heRhoThermo)的构造函数:
    右侧(species::thermo) construct from "dictionary& dict"
    但是mixture model的构造函数例外， i.e. basicMixture(dictionary&, fvMesh&, const word&)
    左侧(basicThermo)
    #1: construct from mesh and phase name
    i.e. rhoThermo(fvMesh&, const word& phaseName);
    #2: construct from mesh, dictionary and phase name
    i.e. rhoThermo(fvMesh&, dictionary&, const word& phaseName); // 这个没用!

//************** pure virtual in basicThermo --- OpenFOAM-2.4.0 **********************//
void correct();
bool incompressible();
bool isochoric();
tmp<volScalarField> rho(); // virtual tmp<scalarField> rho(const label patchi);
volScalarField& he(); // 或const version
tmp<volScalarField> he(const volScalarField& p, const volScalarField& T);
                    // 或he(p, T, const labelList& cells)，或he(p, T, const label patchi)
tmp<volScalarField> hc();
tmp<scalarField> THE(h,p,T0,cells/patchi);
tmp<volScalarField> Cp(); // 或(p, T, patchi)
tmp<volScalarField> Cv(); // 或(p, T, patchi)
tmp<volScalarField> gamma(); // 或(p, T, patchi)
tmp<volScalarField> Cpv(); // 或(p, T, patchi)
tmp<volScalarField> CpByCpv(); // 或(p, T, patchi)

// fields derived from transport state variables
// thermal diffusivity for temperature of mixture (thermal conductivity)
tmp<volScalarField> kappa(); // tmp<scalarField> kappa(const label patchi);

// effective thermal diffusivity for temperature of mixture
tmp<volScalarField> kappaEff(const volScalarField& alphat);
// 或(const scalarField& alphat, const label patchi)

// effective thermal diffusivity of mixture [kg/m/s] 
tmp<volScalarField> alphaEff(alphat); // 或(alphat, patchi)

//**************** pure virtual in fluidThermo --- OpenFOAM-2.4.0 *********************//
volScalarField& psi();
tmp<volScalarField> mu(); // 或(patchi)

//******************* overwrite in psiThermo --- OpenFOAM-2.4.0************************//
tmp<volScalarField> rho(); // virtual tmp<volScalarField> rho(const label patchi);
volScalarField& psi();
tmp<volScalarField> mu(); // 或(patchi)

//******************* overwrite in rhoThermo --- OpenFOAM-2.4.0 ***********************//
tmp<volScalarField> rho(); // virtual tmp<volScalarField> rho(const label patchi);
volScalarField& psi();
volScalarField& mu(); // 或(patchi)

// everithing else overwrite in heThermo!!!

volScalarField&:
p(), he(), T(), alpha() // in basicThermo
psi() // in fluidThermo

tmp<volScalarField>:
// in basicThermo
hc(), Cp(), Cv(), gamma(), Cpv(), CpByCpv(), W(), kappa(), alphahe(), kappaEff(), alphaEff()
mu(), nu() // in fluidThermo



// twoPhaseMixtureThermo in compressibleInterFoam
// 根据压力p和温度T计算内能he, 然后根据压力p和内能he更新每一相的T, psi, rho, mu和alpha
mixture.correctThermo();

// 更新mixture的psi, mu和alpha, 以及更新interfaceProperties (PS: 并没有更新mixture的密度)
// 实际上更新interfaceProperties什么都不会做, 因为表面张力仅与相体积分数有关
// 只是sigma可能与温度有关系
mixture.correct();