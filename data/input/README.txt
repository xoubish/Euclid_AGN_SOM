Astronomy and Astrophysics     Euclid Quick Data Release (Q1): The active Galaxies of Euclid          (Euclid Collaboration: Matamoro Zatarain, T., Fotopoulou, S., Ricci, F., et al., 2025)
================================================================================
    Euclid Quick Data Release (Q1): The active Galaxies of Euclid
    Euclid Collaboration: Matamoro Zatarain, T., Fotopoulou, S., Ricci, F., et al
    <Astron. Astrophys.>
    Accepted
================================================================================
ADC_Keywords: Active gal. nuclei, QSOs, Surveys
Keywords: Galaxies: active, Catalogues, Survey

Abstract: 
We present three catalogues of candidate active galactic nuclei (AGN) in the Euclid Quick Release (Q1) fields. For each Euclid source we collect multi-wavelength photometric and spectroscopic information from surveys such as the Galaxy Evolution Explorer (GALEX), Gaia, Dark Energy Survey (DES), Wide-field Infrared Survey Explorer (WISE), Spitzer, Dark Energy Spectroscopic Instrument (DESI), and Sloan Digital Sky Survey (SDSS), including spectroscopic redshifts from public compilations when available. We investigate the AGN content of the Q1 fields using multiple selection methods. Applying Euclid colours and WISE-AllWISE cuts, we identify 292,222 and 65,131 candidates, respectively. We compile a high-purity QSO catalogue based on Gaia DR3 information, containing 1971 candidates. Using spectroscopic information from DESI, we perform broad-line and narrow-line AGN selections, yielding 4392 AGN candidates across the Q1 fields. We investigate and refine the Euclid Q1 probabilistic random forest QSO population, selecting a refined sample of 180,666 candidates. Additionally, we perform SED fitting on sources with available z_{spec} and, utilizing the derived AGN fraction, identify 7766 AGN candidates. To improve selection purity, we define two new colour criteria (JH_IeY and IeH_gz), finding 313,714 and 267,513 candidates, respectively, across the Q1 fields. We find a total of 229,779 AGN candidates equivalent to an AGN surface density of 3641 deg^{-2} for 18<Ie≤24.5, and a subsample of 30,422 candidates corresponding to an AGN surface density of 482 deg^{-2} when limiting the depth to 18<Ie≤22. The AGN surface densities recovered are consistent with predictions based on AGN X-ray luminosity functions.

Description:
We present three catalogues of Euclid-selected AGN candidates, each corresponding to a separate Euclid field: North, South, and Fornax. The sources have been cross-matched with external surveys including Gaia, GALEX, WISE_AllWISE, DES, SDSS, and DESI. Each entry includes the Euclid source ID, coordinates (RA, Dec), and the corresponding IDs and positions from the matched surveys. When available, spectroscopic redshifts from DESI and SDSS are also included.
In addition to IDs and positional data, the catalogues contain quality flags (good_flags) and classification labels to distinguish stellar and AGN candidates. These labels are based on a range of photometric and spectroscopic selection criteria applied during source identification.

*** NOTE***
To replicate the numbers reported in Euclid Collaboration: Matamoro Zatarain et al. (2025), ensure that QSO/AGN candidates are selected using both the good_flags indicator and the appropriate brightness bin flags (bright_vis_mag_bin, medium_vis_mag_bin, or faint_vis_mag_bin).
************


File Summary:
--------------------------------------------------------------------------------
 FileName     Lrecl  Records  Explanations
--------------------------------------------------------------------------------
ReadMe           80        .  This file
edfn.dat        203        1  Catalogue of sources for the EDF-N
edfs.dat        202        1  Catalogue of sources for the EDF-S
edff.dat        173        1  Catalogue of sources for the EDF-F
--------------------------------------------------------------------------------

	
Byte-by-byte Description of file: edfn.dat
--------------------------------------------------------------------------------
   Bytes Format    Units   Label                            Explanations
--------------------------------------------------------------------------------
   1- 19   I19     ---     object_id_euclid                 Euclid unique source identifier
  21- 37   F17.13  ---     ra_euclid                        Euclid right ascension (J2000.0)
  39- 55   F17.14  ---     dec_euclid                       Euclid declination (J2000.0)
  57- 76   I20     ---     object_id_galex                  GALEX source identifier
  78- 78   E1.0    ---     ra_galex                         GALEX right ascension (J2000.0)
  80- 80   E1.0    ---     dec_galex                        GALEX declination (J2000.0)
  82-101   I20     ---     object_id_gaia                   Gaia unique source identifier
 103-103   E1.0    ---     ra_gaia                          Gaia right ascension (ICRS)
 105-105   E1.0    ---     dec_gaia                         Gaia declination (ICRS)
 107-107   A1      ---     object_id_allwise                WISE-AllWISE unique source identifier
 108-108   E1.0    ---     ra_allwise                       WISE-AllWISE right ascension (J2000.0)
 110-110   E1.0    ---     dec_allwise                      WISE-AllWISE declination (J2000.0)
 112-131   I20     ---     object_id_desi                   DESI unique target ID
 133-133   E1.0    ---     ra_desi                          DESI right ascension (ICRS)
 135-135   E1.0    ---     dec_desi                         DESI declination (ICRS)
 137-137   E1.0    ---     Z_desi                           DESI redshift measured by Redrock
 139-139   A1      ---     object_id_sdss                   SDSS object identification number
 140-140   E1.0    ---     ra_sdss                          SDSS right ascension (J2000.0)
 142-142   E1.0    ---     dec_sdss                         SDSS declination (J2000.0)
 144-144   E1.0    ---     Z_sdss                           SDSS best available redshift taken from Z_VI, Z_PIPE, Z_DR12Q, Z_DR7Q_SCH, Z_DR6Q_HW, and Z_10K
 146-146   I1      ---     good_flags                       Cleaning implemented to keep only those sources with `good flags' (1)
 148-148   I1      ---     bright_vis_mag_bin               Bright Ie magnitude bin: 18<Ie≤21 (1)
 150-150   I1      ---     medium_vis_mag_bin               Medium Ie magnitude bin: 21<Ie≤22 (1)
 152-152   I1      ---     faint_vis_mag_bin                Faint Ie magnitude bin: 22<Ie≤24.5 (1)
 154-154   I1      ---     star_candidate_gaia              Stellar candidate based on Gaia's proper motion and parallax (1)
 156-156   I1      ---     star_candidate_PRF               Stellar candidate based on Probabilistic Random Forest (PRF) star probability > 0.7 (1) 
 158-158   I1      ---     star_candidate_all               Stellar candidate based on star_candidate_gaia and/or star_candidate_PRF (1)
 160-160   I1      ---     PRF_qso_candidate                QSO candidate based on PRF QSO probability > 0.85 (1)
 162-162   I1      ---     B24a_qso_candidate               QSO candidate based on B24A (1)
 164-164   I1      ---     B24b_qso_candidate               QSO candidate based on B24B (1)
 166-166   I1      ---     C75_agn_candidate                AGN candidate based on C75 (1)
 168-168   I1      ---     R90_agn_candidate                AGN candidate based on R90 (1)
 170-170   I1      ---     GDR3_qso_candidate               QSO candidate based on GDR3-QSO (1)
 172-172   I1      ---     JH_IeY_qso_candidate             QSO candidate based on JH_IeY (1)
 174-174   I1      ---     IeH_gz_qso_candidate             QSO candidate based on IeH_gz (1)
 176-176   I1      ---     DESI_qso_candidate               QSO candidate based on DESI SPECTYPE==QSO (1)
 178-178   I1      ---     DESI_broadline_galaxy_candidate  AGN candidate based on DESI SPECTYPE==GALAXY and presence of broad emission lines (1)
 180-180   I1      ---     DESI_broadline_qso_candidate     QSO candidate based on DESI SPECTYPE==QSO and presence of broad emission lines (1)
 182-182   I1      ---     DESI_niibpt_agn_candidate        AGN candidate based on DESI SPECTYPE==GALAXY and Nii BPT diagnostic (1)
 184-184   I1      ---     DESI_siibpt_agn_candidate        AGN candidate based on DESI SPECTYPE==GALAXY and Sii BPT diagnostic (1)
 186-186   I1      ---     DESI_oibpt_agn_candidate         AGN candidate based on DESI SPECTYPE==GALAXY and Oi BPT diagnostic (1)
 188-188   I1      ---     DESI_whan_agn_candidate          AGN candidate based on DESI SPECTYPE==GALAXY and WHAN diagnostic (1)
 190-190   I1      ---     DESI_blue_agn_candidate          AGN candidate based on DESI SPECTYPE==GALAXY and BLUE diagnostic (1)
 192-192   I1      ---     DESI_kex_agn_candidate           AGN candidate based on DESI SPECTYPE==GALAXY and KEX diagnostic (1)
 194-194   E1.0    ---     AGN_fraction                     AGN fraction derived from SED fitting
 196-196   E1.0    ---     AGN_fraction_err                 AGN fraction error derived from SED fitting
 198-203   I6      ---     AGN_sed_candidate                AGN candidate based on AGN fraction > 0.25 (1)
--------------------------------------------------------------------------------
Note (1): True / False values:
	1 = True
	0 = False
--------------------------------------------------------------------------------

Byte-by-byte Description of file: edfs.dat
--------------------------------------------------------------------------------
   Bytes Format    Units   Label                 Explanations
--------------------------------------------------------------------------------
   1- 19   I19     ---     object_id_euclid      Euclid unique source identifier
  21- 36   F16.13  ---     ra_euclid             Euclid right ascension (J2000.0)
  38- 55   F18.14  ---     dec_euclid            Euclid declination (J2000.0)
  57- 76   I20     ---     object_id_galex       GALEX source identifier
  78- 78   E1.0    ---     ra_galex              GALEX right ascension (J2000.0)
  80- 80   E1.0    ---     dec_galex             GALEX declination (J2000.0)
  82-101   I20     ---     object_id_gaia        Gaia unique source identifier
 103-103   E1.0    ---     ra_gaia               Gaia right ascension (ICRS)
 105-105   E1.0    ---     dec_gaia              Gaia declination (ICRS)
 107-126   I20     ---     object_id_des         DES unique identifier for the coadded objects
 128-128   E1.0    ---     ra_des                DES right ascension (J2000.0)
 130-130   E1.0    ---     dec_des               DES declination (J2000.0)
 132-151   A20     ---     object_id_allwise     WISE-AllWISE unique source identifier
 153-162   F10.7   ---     ra_allwise            WISE-AllWISE right ascension (J2000.0)
 164-174   F11.7   ---     dec_allwise           WISE-AllWISE declination (J2000.0)
 176-176   I1      ---     good_flags            Cleaning implemented to keep only those sources with `good flags' (1)
 178-178   I1      ---     bright_vis_mag_bin    Bright Ie magnitude bin: 18<Ie≤21 (1)
 180-180   I1      ---     medium_vis_mag_bin    Medium Ie magnitude bin: 21<Ie≤22 (1)
 182-182   I1      ---     faint_vis_mag_bin     Faint Ie magnitude bin: 22<Ie≤24.5 (1)
 184-184   I1      ---     star_candidate_gaia   Stellar candidate based on Gaia's proper motion and parallax (1)
 186-186   I1      ---     star_candidate_PRF    Stellar candidate based on Probabilistic Random Forest (PRF) star probability > 0.7 (1)
 188-188   I1      ---     star_candidate_all    Stellar candidate based on star_candidate_gaia and/or star_candidate_PRF (1)
 190-190   I1      ---     PRF_qso_candidate     QSO candidate based on PRF QSO probability>0.95 (1)
 192-192   I1      ---     B24a_qso_candidate    QSO candidate based on B24A (1)
 194-194   I1      ---     C75_agn_candidate     AGN candidate based on C75 (1)
 196-196   I1      ---     R90_agn_candidate     AGN candidate based on R90 (1)
 198-198   I1      ---     GDR3_qso_candidate    QSO candidate based on GDR3-QSO (1)
 200-200   I1      ---     JH_IeY_qso_candidate  QSO candidate based on JH_IeY (1)
 202-202   I1      ---     IeH_gz_qso_candidate  QSO candidate based on IeH_gz (1)
--------------------------------------------------------------------------------
Note (1): True / False values:
	1 = True
	0 = False
--------------------------------------------------------------------------------

Byte-by-byte Description of file: edff.dat
--------------------------------------------------------------------------------
   Bytes Format    Units   Label                 Explanations
--------------------------------------------------------------------------------
   1- 19   I19     ---     object_id_euclid      Euclid unique source identifier
  21- 38   F18.15  ---     ra_euclid             Euclid right ascension (J2000.0)
  40- 58   F19.15  ---     dec_euclid            Euclid declination (J2000.0)
  60- 79   I20     ---     object_id_galex       GALEX source identifier
  81- 81   E1.0    ---     ra_galex              GALEX right ascension (J2000.0)
  83- 83   E1.0    ---     dec_galex             GALEX declination (J2000.0)
  85-104   I20     ---     object_id_gaia        Gaia unique source identifier
 106-106   E1.0    ---     ra_gaia               Gaia right ascension (ICRS)
 108-108   E1.0    ---     dec_gaia              Gaia declination (ICRS)
 110-119   I10     ---     object_id_des         DES unique identifier for the coadded objects
 121-129   F9.6    ---     ra_des                DES right ascension (J2000.0)
 131-140   F10.6   ---     dec_des               DES declination (J2000.0)
 142-142   A1      ---     object_id_allwise     WISE-AllWISE unique source identifier
 143-143   E1.0    ---     ra_allwise            WISE-AllWISE right ascension (J2000.0)
 145-145   E1.0    ---     dec_allwise           WISE-AllWISE declination (J2000.0)
 147-147   I1      ---     good_flags            Cleaning implemented to keep only those sources with `good flags' (1)
 149-149   I1      ---     bright_vis_mag_bin    Bright Ie magnitude bin: 18<Ie≤21 (1)
 151-151   I1      ---     medium_vis_mag_bin    Medium Ie magnitude bin: 21<Ie≤22 (1)
 153-153   I1      ---     faint_vis_mag_bin     Faint Ie magnitude bin: 22<Ie≤24.5 (1)
 155-155   I1      ---     star_candidate_gaia   Stellar candidate based on Gaia's proper motion and parallax (1)
 157-157   I1      ---     star_candidate_PRF    Stellar candidate based on Probabilistic Random Forest (PRF) star probability > 0.7 (1)
 159-159   I1      ---     star_candidate_all    Stellar candidate based on star_candidate_gaia and/or star_candidate_PRF (1)
 161-161   I1      ---     PRF_qso_candidate     QSO candidate based on PRF QSO probability>0.95 (1)
 163-163   I1      ---     B24a_qso_candidate    QSO candidate based on B24A (1)
 165-165   I1      ---     C75_agn_candidate     AGN candidate based on C75 (1)
 167-167   I1      ---     R90_agn_candidate     AGN candidate based on R90 (1)
 169-169   I1      ---     GDR3_qso_candidate    QSO candidate based on GDR3-QSO (1)
 171-171   I1      ---     JH_IeY_qso_candidate  QSO candidate based on JH_IeY (1)
 173-173   I1      ---     IeH_gz_qso_candidate  QSO candidate based on IeH_gz (1)
--------------------------------------------------------------------------------
Note (1): True / False values:
	1 = True
	0 = False
--------------------------------------------------------------------------------

Author's address: Teresa Matamoro Zatarain <teresa.matamorozatarain@bristol.ac.uk>

References:

	Cutri, R. M., Wright, E. L., Conrow, T., et al., 2013wise.rept....1C
		Explanatory Supplement to the AllWISE Data Release Products
	Bianchi, L., Shiao, B., & Thilker, D., 2017ApJS..230...24B
		Revised Catalog of GALEX Ultraviolet Sources. I. The All-Sky Survey: GUVcat_AIS
	Assef, R. J., Stern, D., Noirot, G., et al., 2018ApJS..234...23A
		The WISE AGN Catalog
	Abbott, T. M. C., Adamów, M., Aguena, M., et al., 2021ApJS..255...20A
		The Dark Energy Survey Data Release 2
	Abdurro’uf, Accetta, K., Aerts, C., et al. 2022ApJS..259...35A
		The Seventeenth Data Release of the Sloan Digital Sky Surveys: Complete Release of MaNGA, MaStar, and APOGEE-2 Data
	Gaia Collaboration: Vallenari, A., Brown, A. G. A., Prusti, T., et al., 2023A&A...674A...1G
		Gaia Data Release 3. Summary of the content and survey properties
	Euclid Collaboration: Bisigello, L., Massimo, M., Tortora, C., et al., 2024A&A...691A...1E
		Euclid preparation: XLIX. Selecting active galactic nuclei using observed colours
	DESI Collaboration, Adame, A. G., Aguilar, J., et al. 2024AJ....168...58D
		The Early Data Release of the Dark Energy Spectroscopic Instrument
	Storey-Fisher, K., Hogg, D. W., Rix, H.-W., et al., 2024ApJ...964...69S
		Quaia, the Gaia-unWISE Quasar Catalog: An All-sky Spectroscopic Quasar Sample
	Fu, Y., Wu, X.-B., Li, Y., et al. 2024ApJS..271...54F
		CatNorth: An Improved Gaia DR3 Quasar Candidate Catalog with Pan-STARRS1 and CatWISE
	Euclid Collaboration: Tucci, M., Paltani, S., Hartley, W., et al. 2025,A&A submitted, arXiv:2503.15306
		Euclid Quick Data Release (Q1). Photometric redshifts and physical properties of galaxies through the PHZ processing function
	Fu, Y., Wu, X.-B., Bouwens, R. J., et al., arXiv:2503.14141
		The CatSouth Quasar Candidate Catalog for the Southern Sky and a Unified All-Sky Catalog Based on Gaia DR3
	Euclid Collaboration: Aussel, H., Tereno, I., Schirmer, M., et al. 2025, A&A submitted, arXiv:2503.15302
		Euclid Quick Data Release (Q1) -- Data release overview
	
================================================================================
(End)                            Teresa Matamoro Zatarain [University of Bristol, UK]     09-July-2025
