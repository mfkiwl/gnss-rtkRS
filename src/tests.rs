use hifitime::{Epoch, TimeScale};
use nalgebra::Vector3;
use std::str::FromStr;

use crate::{
    orbit::{Keplerian, Orbit, Perturbations},
    prelude::AprioriPosition,
};

use gnss::prelude::{Constellation, SV};

#[test]
fn keplerian_orbit() {
    let sv = SV {
        constellation: Constellation::GPS,
        prn: 1,
    };
    for (epoch, week, secs, kepler, perturb, apriori, elev, azim, ecef) in [(
        Epoch::from_str("2022-01-01T00:00:00 UTC").unwrap(),
        2190,
        518400.0,
        Keplerian {
            a: 26561110.712759566,
            e: 0.00534839148168,
            i_0: 0.957537602313,
            omega_0: 1.03791041521,
            m_0: 2.30316624652,
            omega: -2.3834050415,
        },
        Perturbations {
            dn: 2.3949035344821167e-17,
            i_dot: 5.11807041192e-10,
            omega_dot: -8.0467641439e-09,
            cus: 6.09830021858e-06,
            cuc: 9.85339283943e-07,
            cis: -1.54599547386e-07,
            cic: -1.04308128357e-07,
            crs: 17.3125,
            crc: 258.34375,
        },
        AprioriPosition::from_ecef(Vector3::<f64>::new(3628427.9118, 562059.0936, 5197872.215)),
        8.386332281745226,
        133.44087594021298,
        (16685968.411769923, 20728763.631397538, -1574846.006229475),
    )] {
        let orbit = Orbit::kepler(sv, epoch, week, secs, kepler, perturb, apriori);

        assert_eq!(orbit.sv, sv);

        let err = (
            (ecef.0 - orbit.position.0).abs(),
            (ecef.1 - orbit.position.1).abs(),
            (ecef.2 - orbit.position.2).abs(),
        );

        assert!(err.0 < 1.0E-6, "x(ecef) error too large {}", err.0);
        assert!(err.1 < 1.0E-6, "y(ecef) error too large {}", err.1);
        assert!(err.2 < 1.0E-6, "z(ecef) error too large {}", err.2);

        let elev_err = (orbit.elevation - elev).abs();
        let azim_err = (orbit.azimuth - azim).abs();

        assert!(elev_err < 1.0E-6, "sv_elev err too large {}", elev_err);
        assert!(azim_err < 1.0E-6, "sv_azim error too large {}", azim_err);

        //todo: test other API
        //let orbit = Orbit::position(position);
    }
}
