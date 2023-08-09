//-----------------------------------------------------------------------------//
// cod4X5Y.R                                                                   //
// Logistic model for SWNS grey seal pup production                            //
//                                                                             //
// Copyright 2021 by Landmark Fisheries Research, Ltd.                         //
//                                                                             //
// This software is provided to DFO in the hope that it will be                //
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of      //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                        //
//                                                                             //
// ALL INTELLECTUAL PROPERTY REMAINS WITH LANDMARK FISHERIES RESEARCH, LTD.    //
// THIS SOFTWARE MAY NOT BE REDISTRIBUTED, SUBLICENCED, COPIED, OR SHARED      //
// OUTSIDE OF DFO TECHNOLOGIES WITHOUT THE EXPRESS WRITTEN CONSENT OF          //
// LANDMARK FISHERIES RESEARCH, LTD.                                           //
//                                                                             //
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" //
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   //
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  //
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE    //
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         //
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        //
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    //
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     //
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     //
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  //
// POSSIBILITY OF SUCH DAMAGE.                                                 //
//-----------------------------------------------------------------------------//

#include <TMB.hpp>
#include <iostream>

// isNA - returns TRUE if x is NA and FALSE otherwise
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Objective function ------------------------------------------------------ //
template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(obsI_t);
  DATA_VECTOR(sdI_t);
  DATA_VECTOR(rprior);
  DATA_VECTOR(Kprior);
  
  // Parameters
  PARAMETER(I0);
  PARAMETER(r);
  PARAMETER(K);

  // Dimension sizes
  int nT = obsI_t.size();

  vector<Type> I_t(nT);

  Type objFun = -dnorm( log(r), rprior(0), rprior(1), TRUE );
  objFun -= dnorm( log(K), Kprior(0), Kprior(1), TRUE );

  for( int t=0; t<nT; t++)
  {
    I_t(t) = K*I0/( I0 + (K-I0)*exp(-r*t) );

    if( !isNA(obsI_t(t)) )
      objFun -= dnorm( log(I_t(t)), log(obsI_t(t)), sdI_t(t), TRUE );

  }

  REPORT(I_t);

  return objFun;

}
