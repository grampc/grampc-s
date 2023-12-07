#ifndef VEHICLE_PROBLEM_DESCRIPTION_HPP
#define VEHICLE_PROBLEM_DESCRIPTION_HPP

#include "problem_description/stochastic_problem_description.hpp"
#include <vector>

# define PI           3.14159265358979323846

#define POW2(x) ((x)*(x))

//namespace grampc
//{
    class VehicleProblemDescription : public grampc::StochasticProblemDescription
    {
    public:
        VehicleProblemDescription(const std::vector<typeRNum> &pSys,
                                const std::vector<typeRNum> &pCost,
                                const std::vector<typeRNum> &pCon);

        // Specification of the dimensions
        virtual void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT) override;

        // System dynamics and derivatives
        virtual void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;
        virtual void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p) override;
        virtual void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p) override;
        virtual void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p) override;

        // Integral cost and its derivatives
        virtual void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes) override;
        virtual void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes) override;
        virtual void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes) override;
        virtual void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes) override;

        // Terminal cost and its derivatives
        virtual void Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes) override;
        virtual void dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes) override;
        virtual void dVdp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes) override;
        virtual void dVdT(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes) override;

        // Constraints and its derivatives
        virtual void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;
        virtual void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) override;
        virtual void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) override;
        virtual void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) override;


        /** Additional functions required for Taylor-SMPC */

        /** Jacobian df/dx in vector form (column-wise) **/
        virtual void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;

        /** Jacobian df/dp in vector form (column-wise) **/
        virtual void dfdp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;

        /** Hessian d(df/dx)/dx in vector form (column-wise) */
        virtual void dfdxdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;

        /** Jacobian d(df/dx)/dp in vector form (column-wise) */
        virtual void dfdxdp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;

        /** Hessian d(df/dx)/du in vector form (column-wise) */
        virtual void dfdxdu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;

        /** Jacobian d(df/dp)/du in vector form (column-wise) */
        virtual void dfdpdu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;


        /** Jacobian dh/dx in vector form (column-wise) **/
        virtual void dhdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;

        /** Jacobian dh/du in vector form (column-wise) **/
        virtual void dhdu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;


        /** Hessian d(dh/dx)/dx in vector form (column-wise) **/
        virtual void dhdxdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;

        /** Jacobian d(dh/dx)/du in vector form (column-wise) **/
        virtual void dhdxdu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;
        
    private:
        // system parameters
        std::vector<typeRNum> pSys_;

        // cost parameters
        std::vector<typeRNum> pCost_;

        // constraint parameters
        std::vector<typeRNum> pCon_;

        // vertical forces acting on the front and rear tires
        typeRNum F_f_z_;
        typeRNum F_r_z_;

        // slip angles
        typeRNum alpha_f_;
        typeRNum alpha_r_;

        // lateral forces
        typeRNum F_f_y_;
        typeRNum F_r_y_;

        // temporary values
        typeRNum temp1_;
        typeRNum temp2_;
    };
//}

#endif // VEHICLE_PROBLEM_DESCRIPTION_HPP
