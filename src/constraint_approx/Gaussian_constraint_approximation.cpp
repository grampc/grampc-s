/* This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
 *
 * GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
 *
 * Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
 * All rights reserved.
 *
 * GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 */


#include "constraint_approx/Gaussian_constraint_approximation.hpp"

namespace grampc
{
    GaussianConstraintApproximation::GaussianConstraintApproximation(const Vector& probabilities)
        : integral_vec_{
              0.5000000000, 0.5091099066, 0.5179856166, 0.5298034245, 0.5386407727,
              0.5474143686, 0.5561275953, 0.5676288165, 0.5761528183, 0.5845872442,
              0.5956321593, 0.6037523412, 0.6117572097, 0.6222337919, 0.6299402215,
              0.6399902297, 0.6497922717, 0.6570040188, 0.6663721404, 0.6754949861,
              0.6843758093, 0.6930185998, 0.7014278830, 0.7096085637, 0.7175658024,
              0.7253049208, 0.7328313279, 0.7401504639, 0.7489979664, 0.7558710162,
              0.7625542307, 0.7706435628, 0.7769186758, 0.7845245696, 0.7918789587,
              0.7989508465, 0.8044447458, 0.8110951145, 0.8175263762, 0.8237457188,
              0.8308968920, 0.8366497373, 0.8422135252, 0.8475725313, 0.8538038473,
              0.8588341870, 0.8636781969, 0.8692473291, 0.8745510767, 0.8787762365,
              0.8836466721, 0.8882693661, 0.8927052435, 0.8969206458, 0.9016176344,
              0.9054501965, 0.9097095896, 0.9132005927, 0.9170702022, 0.9207537580,
              0.9242614036, 0.9276026182, 0.9307740268, 0.9337976851, 0.9366812144,
              0.9398096267, 0.9424139779, 0.9452386828, 0.9479151666, 0.9504520622,
              0.9528574161, 0.9551387371, 0.9573030403, 0.9593568863, 0.9615407805,
              0.9633799384, 0.9653379611, 0.9669865318, 0.9687437196, 0.9704041707,
              0.9719692292, 0.9734489324, 0.9750002071, 0.9763126353, 0.9775542913,
              0.9788568462, 0.9800781134, 0.9811162671, 0.9822039232, 0.9832268379,
              0.9841890799, 0.9850943864, 0.9860284334, 0.9868263370, 0.9876483981,
              0.9883504708, 0.9890760915, 0.9897542787, 0.9903896519, 0.9909849819,
              0.9915428533, 0.9920645958, 0.9925967378, 0.9930523351, 0.9935168885,
              0.9939149398, 0.9943207034, 0.9946988932, 0.9950514316, 0.9953801022,
              0.9957110725, 0.9959952024, 0.9962601971, 0.9965272632, 0.9967750669,
              0.9969877072, 0.9972019041, 0.9974007344, 0.9975853164, 0.9977566849,
              0.9979272363, 0.9980739696, 0.9982202191, 0.9983459102, 0.9984713674,
              0.9985870174, 0.9986939855, 0.9987926123, 0.9988838527, 0.9989679947,
              0.9990511286, 0.9991225932, 0.9991931303, 0.9992538484, 0.9993137167,
              0.9993687582, 0.9994193645, 0.9994658952, 0.9995113668, 0.9995504242,
              0.9995863476, 0.9996214894, 0.9996536230, 0.9996812527, 0.9997082770,
              0.9997329938, 0.9997556019, 0.9997775117, 0.9997963263, 0.9998135395,
              0.9998302299, 0.9998454220, 0.9998584626, 0.9998711020, 0.9998826097,
              0.9998930875, 0.9999026279, 0.9999117887, 0.9999196493, 0.9999272052,
              0.9999336834, 0.9999399168, 0.9999455550, 0.9999506675, 0.9999552924,
              0.9999594866, 0.9999632811, 0.9999669025, 0.9999700014, 0.9999729566,
              0.9999754879, 0.9999778999, 0.9999800742, 0.9999820342, 0.9999838013,
              0.9999854725, 0.9999868990, 0.9999881853, 0.9999894028, 0.9999904944,
              0.9999914269, 0.9999923094, 0.9999931007, 0.9999938104, 0.9999944768,
              0.9999950445, 0.9999955537, 0.9999960320, 0.9999964588, 0.9999968223,
              0.9999971636, 0.9999974683, 0.9999977402, 0.9999979829, 0.9999982090,
              0.9999984012, 0.9999985803, 0.9999987325, 0.9999988745, 0.9999990000},
          quantile_vec_{
              0.0000000000, 0.0228371343, 0.0450985380, 0.0747757319, 0.0970099981,
              0.1191313893, 0.1411583993, 0.1703404939, 0.1920610769, 0.2136429395,
              0.2420574067, 0.2630716978, 0.2839018737, 0.3113528208, 0.3316950260,
              0.3584326778, 0.3847597039, 0.4043002219, 0.4299174117, 0.4551379076,
              0.4799704835, 0.5044249335, 0.5285116924, 0.5522415471, 0.5756254202,
              0.5986742100, 0.6213986740, 0.6438093466, 0.6713398289, 0.6930821998,
              0.7145428999, 0.7409679484, 0.7618280306, 0.7875655694, 0.8129580982,
              0.8378796357, 0.8576051897, 0.8819390450, 0.9059784711, 0.9297345089,
              0.9577155450, 0.9807814549, 1.0035968995, 1.0260777496, 1.0528880462,
              1.0750962771, 1.0969949349, 1.1228402633, 1.1481713199, 1.1688910574,
              1.1934155809, 1.2173757959, 1.2410439513, 1.2641987249, 1.2908239494,
              1.3132473690, 1.3389688289, 1.3607306824, 1.3856310382, 1.4101598639,
              1.4343332265, 1.4581661067, 1.4815805167, 1.5046867558, 1.5274966034,
              1.5531774516, 1.5753657546, 1.6003424502, 1.6249666425, 1.6492526998,
              1.6732140375, 1.6968632052, 1.7202119627, 1.7432713495, 1.7688528225,
              1.7913314430, 1.8163018750, 1.8382407579, 1.8626426417, 1.8867670958,
              1.9105569471, 1.9340894138, 1.9599675276, 1.9829346513, 2.0056704731,
              2.0306914160, 2.0553648984, 2.0773695350, 2.1015563876, 2.1254830819,
              2.1491584697, 2.1725894957, 2.1980838206, 2.2210536776, 2.2460125554,
              2.2684949191, 2.2930032637, 2.3172256561, 2.3412227414, 2.3650007138,
              2.3885654931, 2.4118730068, 2.4370770415, 2.4599579428, 2.4846956686,
              2.5071727863, 2.5314687676, 2.5555434094, 2.5794025904, 2.6030519357,
              2.6284357826, 2.6516650627, 2.6747003994, 2.6994473913, 2.7239889577,
              2.7464402642, 2.7705487705, 2.7944672190, 2.8182003049, 2.8417525427,
              2.8668703955, 2.8900291767, 2.9147652413, 2.9375489262, 2.9619206904,
              2.9860607243, 3.0100520875, 3.0338211475, 3.0574492324, 3.0808637770,
              3.1057855767, 3.1288676052, 3.1534118723, 3.1761761059, 3.2003586977,
              3.2243740864, 3.2482257164, 3.2719169287, 3.2969909584, 3.3203173659,
              3.3434940348, 3.3680490703, 3.3924309077, 3.4151390146, 3.4391953357,
              3.4630886364, 3.4868221691, 3.5118634015, 3.5352772829, 3.5585406890,
              3.5831003257, 3.6075046099, 3.6303156677, 3.6543881729, 3.6783143268,
              3.7020968970, 3.7257386156, 3.7505847559, 3.7739214727, 3.7984741018,
              3.8215169558, 3.8457865958, 3.8698776317, 3.8938517775, 3.9176531125,
              3.9413423020, 3.9648639974, 3.9895635532, 4.0128220938, 4.0372266621,
              4.0602311902, 4.0843505597, 4.1083377652, 4.1321950573, 4.1559248158,
              4.1807494249, 4.2041894502, 4.2275086091, 4.2519225378, 4.2761993290,
              4.2991407938, 4.3231570519, 4.3470424487, 4.3707988372, 4.3956043455,
              4.4191040343, 4.4424824249, 4.4669022308, 4.4911977332, 4.5142069342,
              4.5382343612, 4.5621423267, 4.5859342457, 4.6096110040, 4.6342654815,
              4.6576947569, 4.6821104249, 4.7052965177, 4.7294817810, 4.7534243031},
            probabilities_(probabilities),
            coefficients_(probabilities.rows())
    {
        typeInt index_upper;
        typeInt index_lower;
        typeInt j;

        // Go through all constraints
        for (typeInt i = 0; i < probabilities_.rows(); ++i)
        {
            // Find upper and lower bound
            for (j = 1; j < integral_vec_.rows(); ++j)
            {
                if (integral_vec_[j] >= probabilities_[i])
                {
                    break;
                }
            }
            index_upper = j;
            index_lower = j - 1;

            // linear interpolation
            coefficients_[i] = quantile_vec_[index_lower] + (probabilities_[i] - integral_vec_[index_lower]) * (quantile_vec_[index_upper] - quantile_vec_[index_lower]) / (integral_vec_[index_upper] - integral_vec_[index_lower]);
        }
    }

    void GaussianConstraintApproximation::setConstraintProbability(const Vector& probabilities)
    {
        typeInt index_upper;
        typeInt index_lower;
        typeInt j;

        probabilities_ = probabilities;

        // Number of constraints
        typeInt numConstraints = probabilities_.rows();

        // Resize coefficients
        coefficients_.resize(numConstraints);

        // Go through all constraints
        for (typeInt i = 0; i < numConstraints; i++)
        {
            // Find upper and lower bound
            for (j = 1; j < integral_vec_.rows(); ++j)
            {
                if (integral_vec_[j] >= probabilities_[i])
                {
                    break;
                }
            }
            index_upper = j;
            index_lower = j - 1;

            // linear interpolation
            coefficients_[i] = quantile_vec_[index_lower] + (probabilities_[i] - integral_vec_[index_lower]) * (quantile_vec_[index_upper] - quantile_vec_[index_lower]) / (integral_vec_[index_upper] - integral_vec_[index_lower]);
        }
    }

    const Vector& GaussianConstraintApproximation::constraintProbability() const
    {
        return probabilities_;
    }

    const Vector&  GaussianConstraintApproximation::tighteningCoefficient() const
    {
        return coefficients_;
    }

    ChanceConstraintApproximationPtr GaussianApprox(const Vector& probabilities)
    {
        return ChanceConstraintApproximationPtr(new GaussianConstraintApproximation(probabilities));
    }
}