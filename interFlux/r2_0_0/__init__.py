#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV
##


#----------------------------------------------------------------------------
from Foam import ref, man


#----------------------------------------------------------------------------
def _createFields( runTime, mesh ):
        
    ref.ext_Info() << "Reading field p_rgh\n" << ref.nl
    p_rgh = man.volScalarField( man.IOobject( ref.word( "p_rgh" ),
                                              ref.fileName( runTime.timeName() ),
                                              mesh,
                                              ref.IOobject.MUST_READ,
                                              ref.IOobject.AUTO_WRITE ),
                                mesh )
    
    ref.ext_Info() << "Reading field alpha1\n" << ref.nl
    man.interfaceProperties # Load the following library to be able to create corresponding BC - "constantAlphaContactAngleFvPatchScalarField"
    alpha1 = man.volScalarField( man.IOobject( ref.word( "alpha1" ),
                                               ref.fileName( runTime.timeName() ),
                                               mesh,
                                               ref.IOobject.MUST_READ,
                                               ref.IOobject.AUTO_WRITE ),
                                 mesh )
    
    ref.ext_Info() << "Reading field U\n" << ref.nl

    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )
                        
    phi = man.createPhi( runTime, mesh, U )
    
    ref.ext_Info() << "Reading transportProperties\n" << ref.nl
    twoPhaseProperties = man.twoPhaseMixture(U, phi)
    
    rho1 = twoPhaseProperties.rho1()
    rho2 = twoPhaseProperties.rho2()
    
    # Need to store rho for ddt(rho, U)
    rho = man.volScalarField( man.IOobject( ref.word( "rho" ),
                                            ref.fileName( runTime.timeName() ),
                                            mesh,
                                            ref.IOobject.READ_IF_PRESENT ),
                               alpha1 * rho1 + ( 1.0 - alpha1 ) * rho2,
                               alpha1.ext_boundaryField().types() )
    rho.oldTime()
    
    # Mass flux
    # Initialisation does not matter because rhoPhi is reset after the
    # alpha1 solution before it is used in the U equation.
    rhoPhi = man.surfaceScalarField( man.IOobject( ref.word( "rho*phi" ),
                                                   ref.fileName( runTime.timeName() ),
                                                   mesh,
                                                   ref.IOobject.NO_READ,
                                                   ref.IOobject.NO_WRITE ),
                                     rho1 * phi )
    
    # Construct interface from alpha1 distribution
    interface = man.interfaceProperties( alpha1, U, twoPhaseProperties )

    # Construct incompressible turbulence model
    turbulence = man.incompressible.turbulenceModel.New( U, phi, twoPhaseProperties ) 
    
    g = man.readGravitationalAcceleration( runTime, mesh)
    
    #dimensionedVector g0(g);

    #Read the data file and initialise the interpolation table
    #interpolationTable<vector> timeSeriesAcceleration( runTime.path()/runTime.caseConstant()/"acceleration.dat" );
    
    ref.ext_Info() << "Calculating field g.h\n" << ref.nl
    gh = man.volScalarField( ref.word( "gh" ), g & man.volVectorField( mesh.C(), man.Deps( mesh ) ) )
    ghf = man.surfaceScalarField( ref.word( "ghf" ), g & man.surfaceVectorField( mesh.Cf(), man.Deps( mesh ) ) )
    
    p = man.volScalarField( man.IOobject( ref.word( "p" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.NO_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            p_rgh + rho * gh )

    pRefCell = 0
    pRefValue = 0.0
    
    pRefCell, pRefValue = ref.setRefCell( p, p_rgh, mesh.solutionDict().subDict( ref.word( "PIMPLE" ) ), pRefCell, pRefValue )
    
    if p_rgh.needReference():
       p += ref.dimensionedScalar( ref.word( "p" ),
                                   p.dimensions(),
                                   pRefValue - ref.getRefCellValue(p, pRefCell) )
       p_rgh << p - rho * gh
       pass

    return p_rgh, p, alpha1, U, phi, rho1, rho2, rho, rhoPhi, twoPhaseProperties, pRefCell, pRefValue, interface, turbulence, g, gh, ghf
    

#--------------------------------------------------------------------------------------
def setDeltaT( runTime, adjustTimeStep, maxCo, CoNum, maxAlphaCo, alphaCoNum, maxDeltaT ):
    if adjustTimeStep:
       maxDeltaTFact = min( maxCo / ( CoNum + ref.SMALL ), maxAlphaCo / ( alphaCoNum + ref.SMALL ) )
       deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
       
       runTime.setDeltaT( min( deltaTFact * runTime.deltaT().value(), maxDeltaT ) )
       ref.ext_Info() << "deltaT = " <<  runTime.deltaT().value() << ref.nl
       pass
    return runTime


#--------------------------------------------------------------------------------------
def alphaCourantNo( runTime, mesh, alpha1, phi ):
    maxAlphaCo = ref.readScalar( runTime.controlDict().lookup( ref.word( "maxAlphaCo" ) ) )
   
    alphaCoNum = 0.0
    meanAlphaCoNum = 0.0

    if mesh.nInternalFaces():
        tmp = ( alpha1 - 0.01 ).pos() * ( 0.99 - alpha1 ).pos()
        field = ref.fvc.surfaceSum( phi.mag() )
        sumPhi =  tmp() * field.internalField()
       
        alphaCoNum = 0.5 * ( sumPhi / mesh.V().field() ).gMax() * runTime.deltaTValue()
        meanAlphaCoNum = 0.5 * ( sumPhi.gSum() / mesh.V().field().gSum() ) * runTime.deltaTValue()
        pass
  
    ref.ext_Info() << "Interface Courant Number mean: " << meanAlphaCoNum << " max: " << alphaCoNum << ref.nl
    
    return maxAlphaCo, alphaCoNum, meanAlphaCoNum


#----------------------------------------------------------------------------------------
def correctPhi( runTime, mesh, phi, p, p_rgh, rho, U, cumulativeContErr, pimple, pRefCell, pRefValue ):
    
    cumulativeContErr = ref.ContinuityErrs( phi, runTime, mesh, cumulativeContErr )

    pcorrTypes = ref.wordList( p_rgh.ext_boundaryField().size(), ref.zeroGradientFvPatchScalarField.typeName )
    
    for i in range( p.ext_boundaryField().size() ):
       if p_rgh.ext_boundaryField()[i].fixesValue():
          pcorrTypes[i] = ref.fixedValueFvPatchScalarField.typeName
          pass
       pass
    
    pcorr = ref.volScalarField( ref.IOobject( ref.word( "pcorr" ),
                                              ref.fileName( runTime.timeName() ),
                                              mesh,
                                              ref.IOobject.NO_READ,
                                              ref.IOobject.NO_WRITE ),
                            mesh(),
                            ref.dimensionedScalar( ref.word( "pcorr" ), p_rgh.dimensions(), 0.0 ),
                            pcorrTypes  )
    
    rAUf = ref.dimensionedScalar( ref.word( "(1|A(U))" ), ref.dimTime / rho.dimensions(), 1.0)
    
    ref.adjustPhi(phi, U, pcorr)
    
    for nonOrth in range( pimple.nNonOrthCorr() + 1 ):
        pcorrEqn = ref.fvm.laplacian( rAUf, pcorr ) == ref.fvc.div( phi )

        pcorrEqn.setReference(pRefCell, pRefValue)
        pcorrEqn.solve()

        if nonOrth == pimple.nNonOrthCorr():
           phi -= pcorrEqn.flux()
           pass
    
    cumulativeContErr = ref.ContinuityErrs( phi, runTime, mesh, cumulativeContErr )
    
    return cumulativeContErr
    

#--------------------------------------------------------------------------------------
def alphaEqn( mesh, phi, alpha1, rhoPhi, rho1, rho2, interface, nAlphaCorr ):
    alphaScheme = ref.word( "div(phi,alpha)" )
    alpharScheme = ref.word( "div(phirb,alpha)" )
    
    phic = ( phi() / mesh.magSf() ).mag() # mixed calculation
    phic << ( interface.cAlpha() * phic ).ext_min( phic.ext_max() )
    phir = phic * interface.nHatf()
    

    for aCorr in range( nAlphaCorr ):
       phiAlpha = ref.fvc.flux( phi, alpha1, alphaScheme ) + ref.fvc.flux( -ref.fvc.flux( -phir, 1.0 - alpha1, alpharScheme ), alpha1, alpharScheme )
       ref.MULES.explicitSolve( alpha1, phi, phiAlpha, 1.0, 0.0 )
       
       rhoPhi << phiAlpha * ( rho1 - rho2 ) + phi * rho2

       pass
    ref.ext_Info() << "Liquid phase volume fraction = " << alpha1.weightedAverage( mesh.V() ).value() \
               << "  Min(alpha1) = " << alpha1.ext_min().value() \
               << "  Max(alpha1) = " << alpha1.ext_max().value() << ref.nl
    pass


#----------------------------------------------------------------------------------------
def alphaEqnSubCycle( runTime, pimple, mesh, phi, alpha1, rho, rhoPhi, rho1, rho2, interface ):

    nAlphaCorr = ref.readLabel( pimple.dict().lookup( ref.word( "nAlphaCorr" ) ) )
    nAlphaSubCycles = ref.readLabel( pimple.dict().lookup( ref.word( "nAlphaSubCycles" ) ) )
    if (nAlphaSubCycles > 1):
       totalDeltaT = runTime.deltaT()
       rhoPhiSum = 0.0 * rhoPhi

       alphaSubCycle = ref.subCycle_volScalarField(alpha1, nAlphaSubCycles)
       for item in alphaSubCycle: 
           alphaEqn( mesh, phi, alpha1, rhoPhi, rho1, rho2, interface, nAlphaCorr )
           rhoPhiSum += ( runTime.deltaT() / totalDeltaT ) * rhoPhi
           pass
       # To make sure that variable in the local scope will be destroyed
       # - during destruction of this variable it performs some important actions
       # - there is a difference between C++ and Python memory management, namely
       # if C++ automatically destroys stack variables when they exit the scope,
       # Python relay its memory management of some "garbage collection" algorithm
       # that do not provide predictable behavior on the "exit of scope"
       del alphaSubCycle
       
       rhoPhi << rhoPhiSum
    else:
       alphaEqn( mesh, phi, alpha1, rhoPhi, rho1, rho2, interface, nAlphaCorr )
       pass
    interface.correct()
    
    rho == alpha1 * rho1 + ( 1.0 - alpha1 ) * rho2

    pass


#--------------------------------------------------------------------------------------
def _UEqn( mesh, alpha1, U, p, p_rgh, ghf, rho, rhoPhi, turbulence, g, twoPhaseProperties, interface, pimple ):

    muEff = man.surfaceScalarField( ref.word( "muEff" ),
                                    man.surfaceScalarField( twoPhaseProperties.muf(), man.Deps( twoPhaseProperties ) ) + 
                                       man.fvc.interpolate( rho * man.volScalarField( turbulence.ext_nut(), man.Deps( turbulence ) ) ) )
    UEqn = man.fvm.ddt( rho, U ) + man.fvm.div( rhoPhi, U ) - man.fvm.laplacian( muEff, U ) - ( man.fvc.grad( U ) & man.fvc.grad( muEff ) )
    
    UEqn.relax()

    if pimple.momentumPredictor():
       ref.solve( UEqn == \
                    ref.fvc.reconstruct( ( ref.fvc.interpolate( interface.sigmaK() ) * ref.fvc.snGrad( alpha1 ) - ghf * ref.fvc.snGrad( rho ) \
                                                                                                 - ref.fvc.snGrad( p_rgh ) ) * mesh.magSf() ) )
       pass
    
    return UEqn


#--------------------------------------------------------------------------------------
def _pEqn( runTime, mesh, UEqn, U, p, p_rgh, gh, ghf, phi, alpha1, rho, g, interface, corr, pimple, pRefCell, pRefValue, cumulativeContErr ):
    rAU = 1.0 / UEqn.A()
     
    rAUf = ref.fvc.interpolate( rAU )
    
    U << rAU * UEqn.H()
    
    phiU = ref.surfaceScalarField( ref.word( "phiU" ),
                                   ( ref.fvc.interpolate( U ) & mesh.Sf() ) + ref.fvc.ddtPhiCorr( rAU, rho, U, phi ) )
                               
    ref.adjustPhi(phiU, U, p)
    
    phi << ( phiU + ( ref.fvc.interpolate( interface.sigmaK() ) * ref.fvc.snGrad( alpha1 ) - ghf * ref.fvc.snGrad( rho ) )*rAUf*mesh.magSf() )

    for nonOrth in range( pimple.nNonOrthCorr() + 1 ):
        p_rghEqn = ref.fvm.laplacian( rAUf, p_rgh ) == ref.fvc.div( phi ) 
        p_rghEqn.setReference( pRefCell, pRefValue )

        p_rghEqn.solve( mesh.solver( p_rgh.select( pimple.finalInnerIter( corr, nonOrth ) ) ) )
        
        if nonOrth == pimple.nNonOrthCorr():
           phi -= p_rghEqn.flux()
           pass
        pass
    
    U += rAU * ref.fvc.reconstruct( ( phi() - phiU ) / rAUf ) # mixed calculations
    U.correctBoundaryConditions()

    cumulativeContErr = ref.ContinuityErrs( phi, runTime, mesh, cumulativeContErr )
    
    p == p_rgh + rho * gh

    if p_rgh.needReference():
       p += ref.dimensionedScalar( ref.word( "p" ),
                                   p.dimensions(),
                                   pRefValue - ref.getRefCellValue(p, pRefCell) )
       p_rgh << p - rho * gh
       pass
    return cumulativeContErr
    
    
#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )

    mesh = man.createMesh( runTime )
    
    pimple = man.pimpleControl( mesh )
        
    cumulativeContErr = ref.initContinuityErrs()
    
    p_rgh, p, alpha1, U, phi, rho1, rho2, rho, rhoPhi, twoPhaseProperties, pRefCell, \
                                    pRefValue, interface, turbulence, g, gh, ghf = _createFields( runTime, mesh )

    adjustTimeStep, maxCo, maxDeltaT = ref.readTimeControls( runTime )
    
    cumulativeContErr = correctPhi( runTime, mesh, phi, p, p_rgh, rho, U, cumulativeContErr, pimple, pRefCell, pRefValue)
    
    CoNum, meanCoNum = ref.CourantNo( mesh, phi, runTime )
    
    runTime = ref.setInitialDeltaT( runTime, adjustTimeStep, maxCo, CoNum )
    
    ref.ext_Info() << "\nStarting time loop\n" << ref.nl
    
    while runTime.run() :
                
        adjustTimeStep, maxCo, maxDeltaT = ref.readTimeControls( runTime )
        CoNum, meanCoNum = ref.CourantNo( mesh, phi, runTime )        
        maxAlphaCo, alphaCoNum, meanAlphaCoNum = alphaCourantNo( runTime, mesh, alpha1, phi )
        runTime = setDeltaT(  runTime, adjustTimeStep, maxCo, CoNum, maxAlphaCo, alphaCoNum, maxDeltaT )
        
        runTime.increment()
        ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl
        
        twoPhaseProperties.correct()
     
        alphaEqnSubCycle( runTime, pimple, mesh, phi, alpha1, rho, rhoPhi, rho1, rho2, interface )
        
        pimple.start()
        while pimple.loop():
            UEqn = _UEqn( mesh, alpha1, U, p, p_rgh, ghf, rho, rhoPhi, turbulence, g, twoPhaseProperties, interface, pimple )

            # --- PISO loop
            for corr in range( pimple.nCorr() ):
                cumulativeContErr = _pEqn( runTime, mesh, UEqn, U, p, p_rgh, gh, ghf, phi, alpha1, rho, g, 
                                           interface, corr, pimple, pRefCell, pRefValue, cumulativeContErr )
                pass
            
            if pimple.turbCorr():
                turbulence.correct()
                pass
            
            pimple.increment()
            pass
            
        runTime.write()

        ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << \
              "  ClockTime = " << runTime.elapsedClockTime() << " s" << ref.nl << ref.nl
        
        pass

    ref.ext_Info() << "End\n" << ref.nl 

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020000" ):
   if __name__ == "__main__" :
      import sys, os
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
   ref.ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam2.0.0 or higher \n "     
   pass


#--------------------------------------------------------------------------------------

