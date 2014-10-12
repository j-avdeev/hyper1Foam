/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    hyper1Foam

Description
    Solves a transport equation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
 
//#include <armadillo>
#include "complex.H"
 
#include <math.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
// // // //         Info<< "psi before = " << psi << nl << endl;
// // // // 	Info<< "phi before = " << phi << nl << endl;
        
	fvScalarMatrix psiEqn
        (
	    fvm::ddt(psi)
	    +
	    fvm::div(phi,psi)
        );
	



// ===============================================================================

// //                 arma::mat A = arma::zeros<arma::mat>(psi.size(),psi.size());
// // 		arma::mat b = arma::zeros<arma::mat>(psi.size(),psi.size());
// // // 		upSum = p.size();
// //                 forAll(psi,i) // forAll(A,i)
// //                 {
// // //            	    A[i][i] = pEqn.diag()[i];
// //             	    A(i,i) = psiEqn.diag()[i];
// // // // 		    Info << "A("<<i<<","<< i << ") ===" << A(i,i) << nl << endl;
// // //            	    b[i]    = pEqn.source()[i];
// // 		    b(i,0) = psiEqn.source()[i];
// // // // 		    Info << "b("<<i<<","<< 0 << ") ===" << b(i,0) << nl << endl;
// // 		  
// // 		}
// //                 
// //                 const lduAddressing& addr = psiEqn.lduAddr();
// //                 const labelList& lowerAddr = addr.lowerAddr();
// //                 const labelList& upperAddr = addr.upperAddr();                
// //                 
// //                 forAll(lowerAddr, i)
// //                 {
// // 		    A(lowerAddr[i],upperAddr[i]) = psiEqn.upper()[i];
// // // // 		    Info << "A(" << lowerAddr[i] << "," << upperAddr[i] << ")=" << psiEqn.upper()[i] << nl << endl;
// // 		    A(upperAddr[i],lowerAddr[i]) = psiEqn.lower()[i];
// // // // 		    Info << "A(" << upperAddr[i] << "," << lowerAddr[i] << ")=" << psiEqn.lower()[i] << nl << endl;
// // //            	    A[lowerAddr[i]][upperAddr[i]] = pEqn.upper()[i];
// // //            	    A[upperAddr[i]][lowerAddr[i]] = pEqn.lower()[i];           	                	    
// // //             	    downSum += pEqn.upper()[i]* pEqn.upper()[i];
// // //            	    downSum += pEqn.lower()[i]*pEqn.lower()[i];    
// //                 }
// //                 
// //                 forAll(psi.boundaryField(),I) // what is it??
// //                 {
// //             	    const fvPatch &ptch=psi.boundaryField()[I].patch();
// // // // // // // 		    Info << "psi.boundaryField()[" << I << "]" << psi.boundaryField()[I];
// //             	    forAll(ptch,J)
// //             	    {
// //            		int w=ptch.faceCells()[J];
// // // // 			Info << "ptch.faceCells() = " << ptch.faceCells() << nl << endl;
// // //            		A[w][w]+=pEqn.internalCoeffs()[I][J];
// // 			A(w,w)+=psiEqn.internalCoeffs()[I][J];
// // // // 			Info << "A(" << w << "," << w << ") = " << A(w,w) << nl << endl;
// // //            		b[w]   +=pEqn.boundaryCoeffs()[I][J];
// // 			b(w,0) +=psiEqn.boundaryCoeffs()[I][J];
// // // // 			Info << "b(" << w << "," << 0 << ") = " << b(w,0) << nl << endl;
// //             	    }
// //                 
// //                 }
// // // // //                 Info << "psi.boundaryField() = " << psi.boundaryField() << nl << endl;
// // // // // 		Info << "=== A(i,j) ===" << nl << endl;
// // // // // 		for (int i=0; i<psi.size(); i++) // forAll(A,i)
// // // // //                 {
// // // // // 		  for (int j=0; j<psi.size(); j++)
// // // // // 		  {
// // // // // 		    Info << A(i,j) << " ";
// // // // // 		  }
// // // // // 		  Info << nl << endl;
// // // // // 		  
// // // // // 		}
		
		
// // // // // 		Info << "=== b(i,j) ===" << nl << endl;
// // // // // 		for (int i=0; i<psi.size(); i++) // forAll(A,i)
// // // // //                 {
// // // // // // 		  for (int j=0; j<psi.size(); j++)
// // // // // // 		  {
// // // // // 		    Info << b(i,0) << " ";
// // // // // // 		  }
// // // // // 		  Info << nl << endl;
// // // // // 		  
// // // // // 		}
// 		arma::inplace_trans(b);
// // // // // 		arma::mat Ai = arma::inv(A);
// 		arma::mat bt = trans(b);
// // // // // 		arma::mat X = inv(A)*b;
// // 		Info << "=== Ai(i,j) ===" << nl << endl;
// // 		for (int i=0; i<psi.size(); i++) // forAll(A,i)
// //                 {
// // 		  for (int j=0; j<psi.size(); j++)
// // 		  {
// // 		    Info << Ai(i,j) << " ";
// // 		  }
// // 		  Info << nl << endl;
// // 		  
// // 		}
// // 		Info << "=== X(i,j) ===" << nl << endl;
// // 		for (int i=0; i<psi.size(); i++) // forAll(A,i)
// //                 {
// // 		  for (int j=0; j<psi.size(); j++)
// // 		  {
// // 		    Info << X(i,j) << " ";
// // 		  }
// // 		  Info << nl << endl;
// // 		  
// // 		}

// // // // // // // // 		double upSum = 0;			//
// // // // // // // //                 double downSum = 0;			//
// // // // // // // //                 double F = 0;
// // // // // // // //                 arma::mat Ar2 = A;
// // // // // // // // 		arma::inplace_trans(A);
// // // // // // // // 		arma::mat G = A*Ar2;
// // // // // // // // 
// // // // // // // //                 arma::mat D = arma::zeros<arma::mat>(psi.size(),psi.size());
// // // // // // // // 		for (int i=0; i<psi.size(); i++) // forAll(A,i)
// // // // // // // //                 {
// // // // // // // // 		  for (int j=0; j<psi.size(); j++)
// // // // // // // // 		  {
// // // // // // // // 		    if (i==j)
// // // // // // // // 		    {
// 		      if (Gsqr(i,j)==0)
// 		      {
// 			Info << "Gsqr(i,j)==0 !!!" << nl << endl;
// 		      }
// 		      else
// 		      {

// // // // // // // // 			D(i,j) = 1/std::pow(G(i,j),0.5);
// 		      }
// // // // // // // // 		    }
// // // // // // // // 		  }
// // // // // // // // 		}
// // // // // // // // 		arma::mat R = D*G*D;
// // // // // // // // 		upSum = 0;
// // // // // // // // 		for (int i=0; i<psi.size(); i++) // forAll(A,i)
// // // // // // // //                 {
// // // // // // // // 		  for (int j=0; j<psi.size(); j++)
// // // // // // // // 		  {
// // // // // // // // 		    if (i==j)
// // // // // // // // 		    {
// // // // // // // // 		      upSum += R(i,j);
// 		      Info << "R("<< i << "," << j << ")" << R(i,j) << endl;
// 		      Info << "D("<< i << "," << j << ")" << D(i,j) << endl;
// // // // // // // // 		    }
// // // // // // // //         	    downSum += R(i,j)*R(i,j);
// // 		    Info << "downSum = " << downSum << nl << endl;
// // // // // // // // 		  }
// // // // // // // // 		}
// // 		Info << "=== R(i,j) ===" << nl << endl;
// // 		for (int i=0; i<psi.size(); i++) // forAll(A,i)
// //                 {
// // 		  for (int j=0; j<psi.size(); j++)
// // 		  {
// // 		    Info << R(i,j) << " ";
// // 		  }
// // 		  Info << nl << endl;
// // 		  
// // 		}
// // // // // // // // 		Info << "R(" << psi.size() << "," << psi.size() << ")" << R(psi.size()-1,psi.size()-1) << nl << endl;
// // // // // // // //                 scalar x = arma::det(R);
// // // // // // // // 		Info << "det(A)" <<nl<<x<<nl<<endl;
// // // // // // // // 		arma::Col<std::complex<double> > eig = arma::eig_gen(R);
// // // // // // // // 		scalar minEig = min(eig).real();
// // // // // // // // 		scalar maxEig = max(eig).real();
// // // // // // // // 		scalar sumEig = sum(eig).real();
// // // // // // // // 		
// // // // // // // // 		arma::vec eigval;
// // // // // // // // 		arma::mat eigvec;
// // // // // // // // 		arma::eig_sym(eigval, eigvec, R);
// // // // // // // // 		
// // // // // // // // 		Info << "sumEig = " << sumEig << nl << endl;
// // // // // // // // 		if (min(eig).imag()==0 && max(eig).imag()==0)
// // // // // // // // 		Info << "Imaginary parts of minEig and maxEig are zeros" << nl << endl;
// // // // // // // // 		
// // // // // // // // 		Info << "minEig" << nl<< minEig <<nl<<endl;
// // // // // // // // 		Info << "maxEig" << nl<< maxEig <<nl<<endl;
// // // // // // // // 		Info << "Cond^-1" << nl<< minEig/maxEig <<nl<<endl;
// // // // // // // // 		F = upSum*upSum/downSum;
// // // // // // // // 		Info << "upSum=" << upSum << nl << endl;
// // // // // // // // 		Info << "downSum=" << downSum << nl << endl;
// // // // // // // // 		Info << "F=" << F << nl << endl;
// // // // // // // // 		Info << "A(0,0)" << A(0,0) << nl << endl;
		
// ===================================================================================
	
        psiEqn.solve();	
	
// // // // // 	Info<< "psiEqn = " << psiEqn << nl << endl;
// // // // // 	Info<< "psi after = " << psi << nl << endl;
// // // // // 	Info<< "phi after = " << phi << nl << endl;
// // // // // 	Info<< "fvm::ddt(psi) = " << fvm::ddt(psi) << nl << endl;
// // // // // 	Info<< "fvm::div(phi,psi) = " << fvm::div(phi,psi) << nl << endl;
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //