Info<< "Reading field p\n" << endl;

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Sw\n" << endl;

volScalarField Sw
(
    IOobject
    (
        "Sw",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading water dynamic viscosity mu_w\n" << endl;

dimensionedScalar mu_w
(
    transportProperties.lookup("mu_w")
);

Info<< "Reading oil dynamic viscosity mu_o\n" << endl;

dimensionedScalar mu_o
(
    transportProperties.lookup("mu_o")
);


Info<< "Reading porosity prst\n" << endl;

dimensionedScalar prst
(
    transportProperties.lookup("prst")
);


Info<< "Reading absolute permeability K\n" << endl;

dimensionedScalar K
(
    transportProperties.lookup("K")
);

volScalarField krw
(
    IOobject
    (
        "krw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

volScalarField kro
(
    IOobject
    (
        "krw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);


// volScalarField fw
// (
//     IOobject
//     (
//         "fw",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     0.0
// );

// Evaluating alpha value
dimensionedScalar myScalar("myScalar", dimensionSet(0,2,-1,0,0,0,0), 1.0);
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    myScalar
);

volVectorField ut
(
    IOobject
    (
        "ut",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    -(alpha)*fvc::grad(p)
);
ut.write();

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(ut)
);
phi.write();

volScalarField fw
(
    IOobject
    (
        "fw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (Sw/Sw) // Sw
);
fw.write();

#include "createFvModels.H"
// #include "createFvConstraints.H"