// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! This module contains the implementation of dynamic local membrane normal calculation.

use groan_rs::prelude::Vector3D;
use nalgebra::{DMatrix, SVD};

use crate::errors::DynamicNormalError;

/// Calculate membrane normal from a cloud of points using PCA.
///
/// The input is a list of `Vector3D` representing the positions of lipid head groups.
/// Returns an error if fewer than 3 points are provided.
fn membrane_normal_from_cloud(points: &[Vector3D]) -> Result<Vector3D, DynamicNormalError> {
    let num_points = points.len();
    // at least three points are required
    if num_points < 3 {
        return Err(DynamicNormalError::NotEnoughPoints(num_points));
    }

    // calculate the center of the points
    let centroid = points
        .iter()
        .fold(Vector3D::new(0.0, 0.0, 0.0), |sum, p| sum + p)
        / (num_points as f32);

    // demean the coordinates
    let mut data = DMatrix::<f32>::zeros(num_points, 3);
    for (i, p) in points.iter().enumerate() {
        data[(i, 0)] = p.x - centroid.x;
        data[(i, 1)] = p.y - centroid.y;
        data[(i, 2)] = p.z - centroid.z;
    }

    // perform singular value decomposition
    let svd = SVD::new(data, true, true);

    let v_t = match svd.v_t.as_ref() {
        Some(x) => x,
        None => return Err(DynamicNormalError::SVDFailed),
    };
    let n_cols = v_t.nrows();

    // get the last principal direction, which is the normal estimate
    Ok(Vector3D::new(
        v_t[(n_cols - 1, 0)],
        v_t[(n_cols - 1, 1)],
        v_t[(n_cols - 1, 2)],
    )
    .to_unit())
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use groan_rs::{
        prelude::{AtomIteratorWithBox, Sphere},
        system::System,
    };

    use super::*;

    #[test]
    fn test_fail_not_enough_samples() {
        let points = [Vector3D::new(1.0, 2.0, 3.0), Vector3D::new(3.0, 2.0, 1.0)];
        match membrane_normal_from_cloud(&points) {
            Ok(_) => panic!("Function should have failed."),
            Err(DynamicNormalError::NotEnoughPoints(x)) => assert_eq!(x, 2),
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_regular() {
        let mut points = Vec::new();

        // z-axis normal
        for x in [1.0, 2.0, 3.0] {
            for y in [1.0, 2.0, 3.0] {
                points.push(Vector3D::new(x, y, 1.0));
            }
        }

        let normal = membrane_normal_from_cloud(&points).unwrap();
        assert_relative_eq!(normal.x, 0.0);
        assert_relative_eq!(normal.y, 0.0);
        assert_relative_eq!(normal.z.abs(), 1.0);

        points.clear();

        // x-axis normal
        for y in [1.0, 2.0, 3.0] {
            for z in [1.0, 2.0, 3.0] {
                points.push(Vector3D::new(2.0, y, z));
            }
        }

        let normal = membrane_normal_from_cloud(&points).unwrap();
        assert_relative_eq!(normal.x.abs(), 1.0);
        assert_relative_eq!(normal.y, 0.0);
        assert_relative_eq!(normal.z, 0.0);

        points.clear();

        // y-axis normal
        for x in [1.0, 2.0, 3.0] {
            for z in [1.0, 2.0, 3.0] {
                points.push(Vector3D::new(x, -1.5, z));
            }
        }

        let normal = membrane_normal_from_cloud(&points).unwrap();
        assert_relative_eq!(normal.x, 0.0);
        assert_relative_eq!(normal.y.abs(), 1.0);
        assert_relative_eq!(normal.z, 0.0);

        points.clear();

        // tilted
        for xy in [3.0, 4.0, 6.0] {
            for z in [-2.0, -1.0, 0.0] {
                points.push(Vector3D::new(xy, xy, z));
            }
        }

        let normal = membrane_normal_from_cloud(&points).unwrap();
        assert_relative_eq!(normal.x.abs(), 0.7071068);
        assert_relative_eq!(normal.y.abs(), 0.7071068);
        assert_relative_eq!(normal.z, 0.0);
        assert_relative_eq!(normal.x + normal.y, 0.0);
    }

    #[test]
    fn test_slightly_irregular() {
        let mut points = Vec::new();
        let mut zs = [0.96, 0.98, 1.02, 1.03, 1.04, 1.06, 0.99, 0.98, 1.0].into_iter();
        for x in [0.9, 2.1, 3.4] {
            for y in [0.3, 2.2, 3.7] {
                points.push(Vector3D::new(x, y, zs.next().unwrap()));
            }
        }

        let normal = membrane_normal_from_cloud(&points).unwrap();
        assert_relative_eq!(normal.x, 0.0, epsilon = 1e-2);
        assert_relative_eq!(normal.y, 0.0, epsilon = 1e-2);
        assert_relative_eq!(normal.z.abs(), 1.0, epsilon = 1e-2);
    }

    #[test]
    fn test_real_planar() {
        let mut system = System::from_file("tests/files/pcpepg.tpr").unwrap();
        system.group_create("Phosphori", "name P").unwrap();

        let expected_normals = [
            Vector3D::new(0.019939683, -0.07379772, -0.9970739),
            Vector3D::new(0.03537177, 0.1589833, 0.9866475),
            Vector3D::new(0.1259833, 0.11300544, 0.98557496),
            Vector3D::new(0.10208496, 0.037488986, 0.99406904),
            Vector3D::new(-0.1796861, -0.13556388, -0.9743385),
            Vector3D::new(0.08554515, -0.063694805, -0.99429625),
            Vector3D::new(0.058934886, -0.048600454, -0.9970781),
            Vector3D::new(-0.15906493, -0.059976142, 0.9854447),
            Vector3D::new(0.07179594, -0.07779022, 0.99438125),
            Vector3D::new(-0.16887201, -0.03536218, 0.9850035),
            Vector3D::new(0.054649983, 0.17095369, -0.9837623),
            Vector3D::new(0.06964987, -0.068960726, -0.9951851),
            Vector3D::new(0.12355773, 0.0035226601, -0.99233115),
            Vector3D::new(0.14382878, -0.03739055, -0.988896),
            Vector3D::new(0.08394643, 0.0697191, -0.99402833),
            Vector3D::new(0.056339037, 0.0075176265, -0.9983834),
            Vector3D::new(0.07272931, -0.014148992, -0.99725133),
            Vector3D::new(-0.11291739, 0.06325698, -0.9915888),
            Vector3D::new(0.043125153, 0.008952808, 0.9990296),
            Vector3D::new(-0.16973464, 0.054992236, -0.98395425),
            Vector3D::new(0.030309448, 0.03769139, 0.99882966),
            Vector3D::new(0.05035916, -0.0010272319, 0.9987307),
            Vector3D::new(0.019779278, -0.008824456, -0.99976546),
            Vector3D::new(0.02662054, 0.0029694305, 0.99964124),
            Vector3D::new(0.10956673, -0.030122284, -0.99352294),
            Vector3D::new(0.14183103, -0.05797418, 0.98819184),
            Vector3D::new(0.03373865, 0.04190623, 0.9985518),
            Vector3D::new(-0.14346373, 0.06971726, -0.9871969),
            Vector3D::new(0.027205076, 0.061963305, 0.9977076),
            Vector3D::new(0.017223975, 0.07523002, 0.9970175),
            Vector3D::new(0.018885745, 0.01152975, -0.9997552),
            Vector3D::new(0.002023952, -0.014715449, 0.99988973),
            Vector3D::new(0.1099667, 0.025607673, -0.9936054),
            Vector3D::new(0.007469739, -0.026535803, 0.99961996),
            Vector3D::new(0.011718081, 0.036759418, 0.9992554),
            Vector3D::new(0.05790448, 0.025840282, -0.9979877),
            Vector3D::new(-0.18638445, 0.072998166, -0.97976124),
            Vector3D::new(0.051643915, 0.01804022, -0.9985027),
            Vector3D::new(0.03453225, -0.007023326, 0.9993789),
            Vector3D::new(0.017506903, -0.04317353, 0.9989142),
            Vector3D::new(0.018763635, -0.01350592, 0.99973273),
            Vector3D::new(0.016533308, 0.023119433, 0.999596),
            Vector3D::new(-0.1356995, 0.08210372, -0.98734224),
            Vector3D::new(0.0662576, 0.003040572, -0.99779797),
            Vector3D::new(0.12003526, 0.0039590574, 0.99276173),
            Vector3D::new(0.03380418, 0.0125501165, -0.9993497),
            Vector3D::new(0.17817241, -0.023221834, 0.98372525),
            Vector3D::new(0.014899884, 0.06330999, 0.9978827),
            Vector3D::new(0.008362867, 0.005710643, -0.99994874),
            Vector3D::new(0.0125309955, -0.010883104, 0.99986225),
            Vector3D::new(0.007830463, 0.026867194, -0.9996084),
            Vector3D::new(0.076697186, 0.023985201, -0.9967659),
            Vector3D::new(0.020920651, 0.012544673, 0.9997025),
            Vector3D::new(0.017652782, -0.005117607, 0.9998311),
            Vector3D::new(0.07683702, 0.027184965, -0.99667305),
            Vector3D::new(-0.019220488, 0.014189856, -0.9997146),
            Vector3D::new(0.07167118, -0.0076096975, 0.99739933),
            Vector3D::new(0.102221645, 0.01109694, 0.9946998),
            Vector3D::new(0.021864213, -0.005639493, 0.99974513),
            Vector3D::new(0.013446753, -0.004123769, 0.9999011),
            Vector3D::new(0.020246007, -0.01450603, 0.9996898),
            Vector3D::new(0.104586385, 0.0013575138, 0.9945149),
            Vector3D::new(0.005310551, 0.016070712, -0.99985677),
            Vector3D::new(0.018742599, -0.009936025, 0.999775),
            Vector3D::new(0.011090281, -0.001454881, 0.9999375),
            Vector3D::new(-0.1171862, -0.058579367, -0.9913809),
            Vector3D::new(0.1377209, 0.051987126, 0.9891058),
            Vector3D::new(0.13688965, 0.056458835, 0.9889761),
            Vector3D::new(0.09769562, 0.053360034, 0.9937849),
            Vector3D::new(0.15269312, -0.021383153, 0.9880423),
            Vector3D::new(-0.045325033, -0.03272814, -0.9984361),
            Vector3D::new(0.22422183, 0.06659179, 0.9722603),
            Vector3D::new(0.0656786, 0.11274521, 0.9914509),
            Vector3D::new(0.028027311, 0.0020895798, 0.99960494),
            Vector3D::new(0.012348495, -0.0019297748, 0.9999219),
            Vector3D::new(0.18561861, 0.049006406, 0.98139906),
            Vector3D::new(0.17894161, 0.07470498, 0.98101944),
            Vector3D::new(0.20964426, 0.020942084, -0.9775534),
            Vector3D::new(0.01450741, -0.045000173, 0.9988817),
            Vector3D::new(0.029621119, 0.069938645, -0.99711144),
            Vector3D::new(0.01737121, 0.0129034575, -0.9997659),
            Vector3D::new(0.05136815, 0.014987842, 0.99856734),
            Vector3D::new(0.058922336, -0.03741177, 0.99756134),
            Vector3D::new(0.03059999, -0.09067008, 0.9954108),
            Vector3D::new(0.014156878, 0.03503776, 0.99928576),
            Vector3D::new(0.083881535, -0.0072715236, 0.99644923),
            Vector3D::new(0.010266017, 0.09858759, -0.99507546),
            Vector3D::new(0.06673409, 0.0050728694, 0.997758),
            Vector3D::new(0.1316364, 0.03357398, 0.9907294),
            Vector3D::new(0.007314736, 0.009987437, 0.9999234),
            Vector3D::new(0.05648798, 0.06879973, -0.99603003),
            Vector3D::new(0.0115753375, -0.054214403, 0.9984622),
            Vector3D::new(0.025158277, 0.009130374, 0.9996418),
            Vector3D::new(0.027361248, 0.010647625, 0.99956894),
            Vector3D::new(0.18030128, -0.040465247, -0.9827788),
            Vector3D::new(0.060588676, 0.014990018, 0.99805033),
            Vector3D::new(0.001503393, 0.0022171568, -0.9999964),
            Vector3D::new(0.10523872, 0.010429123, 0.9943923),
            Vector3D::new(0.050380222, 0.071050934, 0.9961996),
            Vector3D::new(0.08738596, 0.053290714, -0.9947481),
            Vector3D::new(0.0062247305, 0.022175215, 0.99973476),
            Vector3D::new(0.0071625393, -0.044201687, 0.99899703),
            Vector3D::new(0.064999625, -0.01144854, 0.9978196),
            Vector3D::new(0.012825869, 0.04117396, 0.9990697),
            Vector3D::new(0.00856402, -0.057089534, 0.9983323),
            Vector3D::new(0.06869478, 0.020969523, 0.99741733),
            Vector3D::new(0.10237356, 0.020989606, 0.9945246),
            Vector3D::new(0.10747415, 0.0057116873, -0.99419147),
            Vector3D::new(0.0058913156, 0.022786628, -0.999723),
            Vector3D::new(0.005185337, -0.080362804, 0.9967522),
            Vector3D::new(0.011266603, -0.017862618, -0.999777),
            Vector3D::new(0.0029963565, 0.014647834, 0.99988824),
            Vector3D::new(0.04381818, 0.022273785, 0.9987912),
            Vector3D::new(0.038245603, 0.018197304, -0.99910265),
            Vector3D::new(0.0032856727, 0.014930252, 0.9998832),
            Vector3D::new(0.02165919, 0.0018291819, 0.9997637),
            Vector3D::new(0.012469897, -0.0024654237, -0.99991924),
            Vector3D::new(0.07369235, 0.024658548, 0.99697614),
            Vector3D::new(0.026197908, 0.0013273709, 0.9996559),
            Vector3D::new(0.022671504, -0.0016797365, 0.99974155),
            Vector3D::new(0.032641158, 0.007563734, 0.9994385),
            Vector3D::new(0.005198498, -0.012342119, 0.99991035),
            Vector3D::new(0.014842956, -0.103376724, -0.9945316),
            Vector3D::new(0.013717366, 0.06675895, 0.9976748),
            Vector3D::new(0.01934048, 0.0031652905, 0.99980795),
            Vector3D::new(0.025475524, -0.0146073885, -0.99956876),
            Vector3D::new(0.081691496, 0.014502645, -0.9965522),
            Vector3D::new(0.0064702774, -0.003546236, 0.9999728),
            Vector3D::new(0.013282725, 0.06777206, 0.9976124),
            Vector3D::new(0.08629104, 0.0062431656, -0.99625045),
            Vector3D::new(0.023354126, -0.008517809, 0.99969095),
            Vector3D::new(0.16517706, 0.0451389, -0.98523045),
            Vector3D::new(0.067722976, 0.085232615, -0.9940568),
            Vector3D::new(0.07779883, -0.10074108, -0.99186623),
            Vector3D::new(-0.19473654, -0.08404782, 0.977248),
            Vector3D::new(0.057394642, 0.08744711, 0.9945144),
            Vector3D::new(0.01841816, 0.14645308, -0.98904616),
            Vector3D::new(0.057153635, -0.009662875, -0.9983187),
            Vector3D::new(0.020236697, -0.0739179, -0.99705905),
            Vector3D::new(0.1639796, 0.08627973, 0.98268336),
            Vector3D::new(0.14045754, 0.06749089, -0.98778373),
            Vector3D::new(0.13182275, 0.11170862, -0.9849589),
            Vector3D::new(0.12890683, 0.04877548, -0.99045646),
            Vector3D::new(0.07964923, 0.00836785, 0.99678785),
            Vector3D::new(0.16063029, -0.008960532, 0.986974),
            Vector3D::new(0.05818244, -0.047678437, -0.99716675),
            Vector3D::new(-0.19879247, -0.016502095, -0.9799027),
            Vector3D::new(0.09642512, 0.04351972, -0.9943884),
            Vector3D::new(-0.18673299, -0.034614902, -0.98180073),
            Vector3D::new(0.07631385, 0.0871758, -0.9932656),
            Vector3D::new(0.011463047, -0.060175683, 0.99812204),
            Vector3D::new(0.073577054, -0.13329926, 0.9883409),
            Vector3D::new(0.14249206, -0.11049856, 0.9836087),
            Vector3D::new(0.028860854, -0.02402097, -0.9992948),
            Vector3D::new(-0.12580992, -0.048694313, 0.9908586),
            Vector3D::new(0.07471252, 0.03274562, -0.9966673),
            Vector3D::new(0.1780791, 0.055236466, 0.9824647),
            Vector3D::new(0.06640337, 0.01364299, -0.9976996),
            Vector3D::new(0.07230141, -0.010115322, -0.99733156),
            Vector3D::new(0.008830541, -0.01860479, 0.999788),
            Vector3D::new(0.05277334, 0.013062054, -0.99852115),
            Vector3D::new(0.03703564, 0.0048170164, -0.9993023),
            Vector3D::new(0.04154935, 0.032464318, 0.99860895),
            Vector3D::new(0.0029357865, -0.049153887, 0.9987869),
            Vector3D::new(0.014484387, 0.014592433, 0.99978864),
            Vector3D::new(0.008376108, 0.020997116, -0.9997445),
            Vector3D::new(0.03448549, 0.002914466, -0.999401),
            Vector3D::new(-0.16898082, 0.10744703, -0.97974515),
            Vector3D::new(0.13364373, -0.038943343, -0.99026406),
            Vector3D::new(0.09022468, -0.0000849648, 0.99592143),
            Vector3D::new(0.032458894, 0.039232254, 0.9987028),
            Vector3D::new(0.0015771651, 0.06424951, -0.9979326),
            Vector3D::new(0.01770657, 0.011696967, 0.9997748),
            Vector3D::new(0.06025689, 0.015405527, -0.99806404),
            Vector3D::new(0.0015214644, 0.009084954, 0.99995756),
            Vector3D::new(0.0027014478, -0.05827556, -0.99829686),
            Vector3D::new(0.019031407, -0.0038886494, 0.99981135),
            Vector3D::new(0.052264825, 0.019872546, -0.99843556),
            Vector3D::new(0.0012935586, 0.053675197, -0.9985576),
            Vector3D::new(0.012418912, -0.019796092, 0.99972695),
            Vector3D::new(0.014450457, -0.098967984, 0.9949857),
            Vector3D::new(0.008271486, 0.04539237, 0.998935),
            Vector3D::new(0.000885558, 0.021643898, -0.99976534),
            Vector3D::new(0.0127518475, -0.00498249, 0.99990636),
            Vector3D::new(0.0853087, 0.00089619064, 0.99635416),
            Vector3D::new(0.0077899448, -0.0038013898, 0.99996245),
            Vector3D::new(0.011018865, -0.0110073695, 0.9998787),
            Vector3D::new(0.010936995, 0.008424302, 0.9999047),
            Vector3D::new(0.0029810355, -0.028259093, 0.9995962),
            Vector3D::new(0.13784367, -0.002339523, 0.9904513),
            Vector3D::new(0.011435856, 0.00048372388, 0.99993455),
            Vector3D::new(0.014019656, -0.009798379, 0.99985373),
            Vector3D::new(0.016148085, -0.08076457, 0.9966024),
            Vector3D::new(0.010630005, 0.0051639746, 0.9999302),
            Vector3D::new(0.021436807, -0.0017171323, 0.99976873),
            Vector3D::new(0.124009274, 0.032204166, 0.9917584),
            Vector3D::new(0.07402569, 0.028147593, 0.996859),
            Vector3D::new(0.075847164, 0.025947856, 0.99678177),
            Vector3D::new(0.08849525, 0.11413847, 0.98951554),
            Vector3D::new(0.054100893, 0.020795844, 0.9983189),
            Vector3D::new(0.031789456, 0.0397228, 0.99870497),
            Vector3D::new(0.11645165, 0.05523297, 0.9916594),
            Vector3D::new(0.04645506, -0.037698504, 0.99820876),
            Vector3D::new(0.20632994, 0.037386943, 0.977768),
            Vector3D::new(0.02490566, -0.004060921, 0.99968153),
            Vector3D::new(-0.1937008, -0.067636825, -0.9787264),
            Vector3D::new(0.0053090234, -0.009488632, 0.99994093),
            Vector3D::new(0.04165643, -0.009548769, 0.9990864),
            Vector3D::new(0.119337976, 0.022006009, 0.99260986),
            Vector3D::new(0.07294969, -0.035910457, 0.99668896),
            Vector3D::new(0.15078378, 0.08898306, 0.9845538),
            Vector3D::new(0.044224415, -0.05671286, 0.9974106),
            Vector3D::new(0.04663332, 0.01975616, 0.9987167),
            Vector3D::new(0.17334202, 0.058827203, 0.9831032),
            Vector3D::new(0.008287059, -0.028548414, 0.9995581),
            Vector3D::new(0.15930898, 0.03247035, 0.9866947),
            Vector3D::new(0.001636836, 0.005889465, 0.99998134),
            Vector3D::new(0.048841543, 0.026911909, 0.9984439),
            Vector3D::new(0.048187852, -0.06501943, -0.99671984),
            Vector3D::new(0.043974694, 0.02533801, -0.99871135),
            Vector3D::new(0.007472532, 0.0028095604, 0.9999682),
            Vector3D::new(-0.1119736, -0.09400193, -0.9892551),
            Vector3D::new(0.0029338957, 0.012881504, 0.99991274),
            Vector3D::new(0.13301614, 0.0077046542, 0.9910839),
            Vector3D::new(0.021625997, 0.022194874, 0.99951977),
            Vector3D::new(0.00013178881, -0.05760292, 0.9983396),
            Vector3D::new(0.0026647332, 0.058530964, 0.998282),
            Vector3D::new(0.010255866, 0.059680574, -0.99816483),
            Vector3D::new(0.045624543, -0.015945725, 0.9988314),
            Vector3D::new(0.017341733, 0.052962717, 0.9984459),
            Vector3D::new(0.10209441, 0.009004121, 0.994734),
            Vector3D::new(0.018687421, 0.076155744, 0.9969208),
            Vector3D::new(0.00008859144, -0.031789158, 0.9994946),
            Vector3D::new(0.11221505, 0.0113253165, 0.99361944),
            Vector3D::new(0.0012944645, 0.019485723, 0.9998093),
            Vector3D::new(0.012737949, -0.051242243, 0.998605),
            Vector3D::new(0.00087471365, -0.06773704, 0.99770284),
            Vector3D::new(0.13040422, 0.008042491, 0.9914283),
            Vector3D::new(0.011476024, 0.05274674, 0.998542),
            Vector3D::new(-0.18363026, 0.026091792, 0.982649),
            Vector3D::new(0.002241636, 0.046172317, 0.998931),
            Vector3D::new(0.01111322, 0.013121705, -0.9998522),
            Vector3D::new(0.008486928, -0.024959296, 0.99965245),
            Vector3D::new(0.0070460327, 0.020530747, 0.9997644),
            Vector3D::new(0.02918306, -0.0014788262, 0.999573),
            Vector3D::new(0.04034323, -0.003967859, -0.99917805),
            Vector3D::new(0.00077294686, 0.0051418347, -0.99998647),
            Vector3D::new(0.04794095, -0.015608487, -0.9987282),
            Vector3D::new(0.00083597604, 0.01297366, -0.9999155),
            Vector3D::new(0.0007716491, 0.016444162, -0.9998645),
            Vector3D::new(0.019869087, 0.0660345, 0.9976195),
            Vector3D::new(0.021286836, 0.0039114333, 0.9997658),
            Vector3D::new(0.020334888, 0.002284708, 0.9997906),
            Vector3D::new(0.00005519746, -0.00796756, 0.9999683),
            Vector3D::new(0.012352425, -0.0046472875, 0.9999129),
            Vector3D::new(0.0015263354, -0.013447059, -0.99990845),
            Vector3D::new(0.017798379, 0.0046646725, 0.9998308),
            Vector3D::new(0.019559924, 0.0021126256, 0.9998065),
            Vector3D::new(0.017250769, -0.0021558127, 0.9998489),
            Vector3D::new(-0.13809898, 0.035284396, -0.9897897),
            Vector3D::new(0.22083761, 0.09694856, 0.97048014),
            Vector3D::new(0.0028000707, 0.008177582, -0.9999626),
            Vector3D::new(0.056654707, 0.0038809413, -0.9983863),
            Vector3D::new(0.012180123, 0.01637695, -0.99979174),
            Vector3D::new(0.07223803, -0.0058639552, 0.9973702),
            Vector3D::new(0.017795023, 0.035873506, 0.9991979),
            Vector3D::new(0.009829432, 0.053457435, 0.99852175),
            Vector3D::new(0.10387018, 0.04247612, -0.9936835),
            Vector3D::new(0.01344039, 0.020788163, -0.9996936),
            Vector3D::new(0.008888947, 0.0020525972, 0.9999584),
            Vector3D::new(0.13179842, 0.010032575, 0.9912258),
            Vector3D::new(0.0005216558, 0.05326568, -0.9985803),
            Vector3D::new(0.05637889, 0.017959261, 0.9982479),
            Vector3D::new(0.009247412, -0.033216994, 0.9994054),
        ];

        for (i, atom) in system.group_iter("Phosphori").unwrap().enumerate() {
            let geometry = Sphere::new(atom.get_position().unwrap().clone(), 2.0);
            let positions = system
                .group_iter("Phosphori")
                .unwrap()
                .filter_geometry(geometry)
                .map(|a| a.get_position().unwrap().clone())
                .collect::<Vec<Vector3D>>();
            let normal = membrane_normal_from_cloud(&positions).unwrap();

            assert_relative_eq!(normal.x, expected_normals[i].x);
            assert_relative_eq!(normal.y, expected_normals[i].y);
            assert_relative_eq!(normal.z, expected_normals[i].z);
        }
    }
}
