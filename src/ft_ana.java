import java.io.*;

import org.jlab.clas.physics.LorentzVector;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;


/**
 * @author akhanal
 *
 */

public class ft_ana {
	//public int NFTElec;
	public float EB, Eb, Mp;
	public float STT, RFT, FTSTT;
	
	
	public LorentzVector VB, VT, Ve, VGS, Vprot, Vpip;
	
	public int e_ft_part_ind;
	
	public float e_mom, e_the, e_phi;
	public float e_xB, e_Q2, e_W;
	
	public int prot_part_ind;
	public float prot_mom, prot_the, prot_phi, prot_vx, prot_vy, prot_vz, prot_beta, prot_ftb_beta;
	
	public int pip_part_ind, pip_FTOF_pad1b;
	public float pip_mom, pip_the, pip_phi, pip_vx, pip_vy, pip_vz, pip_ftb_beta, pip_FTOF1b_t, pip_FTOF1b_path, pip_FTOF1b_vt;
	
	public H1F H_FT_W;
	public H2F H_FT_e_beta_mom;
	public H2F H_FT_e_t_f, H_FT_e_p_f, H_FT_e_p_the;
	public H2F H_FT_W_Q2, H_FT_e_xB_Q2;
	public H2F[] H_FTOF_pos_beta_mom_pad1a, H_FTOF_neg_beta_mom_pad1a, H_FTOF_pos_beta_mom_pad1b, H_FTOF_neg_beta_mom_pad1b;
	public H2F[] H_FTOF_pos_mass_mom_pad1a, H_FTOF_pos_mass_the_pad1a, H_FTOF_neg_mass_mom_pad1a, H_FTOF_neg_mass_the_pad1a;
	public H2F[] H_FTOF_pos_mass_mom_pad1b, H_FTOF_pos_mass_the_pad1b, H_FTOF_neg_mass_mom_pad1b, H_FTOF_neg_mass_the_pad1b;
	public H2F H_FD_pos_beta_mom, H_FD_neg_beta_mom, H_FD_neutral_beta_mom;
	public H2F H_FD_pos_mass_mom, H_FD_neg_mass_mom, H_FD_neutral_mass_mom;
	public H2F H_FD_pos_mass_the, H_FD_neg_mass_the, H_FD_neutral_mass_the;
	public H2F H_FD_pos_mass_phi, H_FD_neg_mass_phi, H_FD_neutral_mass_phi;
	public F1D F_prot_beta_mom, F_kp_beta_mom, F_pip_beta_mom; 
//	public H2F[] H_FD_pos_mass_mom_pad1b, H_FD_pos_mass_the_pad1b, H_FD_neg_mass_mom_pad1b, H_FD_neg_mass_the_pad1b;
	
	
	public ft_ana() {
	//	NFTElec = 0;
		Eb = 10.6f;
		Mp = 0.93827f;
		
		VB = new LorentzVector(0, 0, Eb, Eb);
		VT = new LorentzVector(0, 0, 0, Mp);
		
		// theoretical 1D functions for prot, kaon and pion
		
		F_prot_beta_mom = new F1D("F_prot_beta_mom", "x/sqrt(0.93827*0.93827+x*x)", 0.3, 5.0);
		F_prot_beta_mom.setLineWidth(2);
		F_prot_beta_mom.setLineColor(2);
		F_kp_beta_mom = new F1D("F_kp_beta_mom", "x/sqrt(0.49367*0.49367+x*x)", 0.3, 5.0);
		F_kp_beta_mom.setLineColor(2);
		F_kp_beta_mom.setLineWidth(2);
		F_pip_beta_mom = new F1D("F_pip_beta_mom", "x/sqrt(0.13957*0.13957+x*x)", 0.3, 5.0);
		F_pip_beta_mom.setLineColor(2);
		F_pip_beta_mom.setLineWidth(2);
		
		H_FT_e_t_f = new H2F("H_FT_e_t_f", "H_FT_e_t_f", 100, -180, 180, 50, 0, 8);
		H_FT_e_t_f.setTitle("electron #theta vs #phi");
		H_FT_e_t_f.setTitleX("#phi (^o)");
		H_FT_e_t_f.setTitleY("#theta (^o)");
		
		H_FT_e_p_the = new H2F("H_FT_e_p_the", "H_FT_e_p_the", 100, 0, 8, 100, 0, 12);
		H_FT_e_p_the.setTitle("electron p vs #theta (^o)");
		H_FT_e_p_the.setTitleX("#theta (^o)");
		H_FT_e_p_the.setTitleY("p (GeV)");

		H_FT_e_p_f = new H2F("H_FT_e_p_f", "H_FT_e_p_f", 100, -180, 180, 100, 0, 12);
		H_FT_e_p_f.setTitle("electron p vs #phi");
		H_FT_e_p_f.setTitleX("#phi (^o)");
		H_FT_e_p_f.setTitleY("p (GeV)");
		
		H_FT_W_Q2 = new H2F("H_FT_W_Q2", "H_FT_W_Q2", 100, 0, 5, 100, 0.00001, 0.7);
		H_FT_W_Q2.setTitle("FT Q^2 vs W");
		H_FT_W_Q2.setTitleX("W ( GeV)");
		H_FT_W_Q2.setTitleY("Q^2 (GeV^2)");
		
		H_FT_e_xB_Q2 = new H2F("H_FT_e_xB_Q2", "H_FT_e_xB_Q2", 100, 0, 0.35, 100, 0.00001, 0.7);
		H_FT_e_xB_Q2.setTitle("electron Q^2 vs xB");
		H_FT_e_xB_Q2.setTitleX("xB");
		H_FT_e_xB_Q2.setTitleY("Q^2 (GeV^2)");
		
		H_FT_W = new H1F("H_FT_W", "H_FT_W", 100, 0, 5);
		H_FT_W.setTitle("electron W");
		H_FT_W.setTitleX("W (GeV)");
		H_FT_W.setTitleY("count");
		
		H_FT_e_beta_mom = new H2F("H_FT_e_beta_mom", "H_FT_e_beta_mom", 100, 0, 10.6, 100, 0, 1.2);
		H_FD_pos_beta_mom = new H2F("H_FD_pos_beta_mom", "H_FD_pos_beta_mom", 100, 0, 10.6, 100, 0, 1.2);
		H_FD_pos_beta_mom.setTitle("POS  #beta vs mom");
		H_FD_pos_beta_mom.setTitleX("p (GeV)");
		H_FD_pos_beta_mom.setTitleY("FTB #beta");
		H_FD_neg_beta_mom = new H2F("H_FD_neg_beta_mom", "H_FD_neg_beta_mom", 100, 0, 10.6, 100, 0, 1.2);
		H_FD_neg_beta_mom.setTitle("NEG  #beta vs mom");
		H_FD_neg_beta_mom.setTitleX("p (GeV)");
		H_FD_neg_beta_mom.setTitleY("FTB #beta");
		H_FD_neutral_beta_mom = new H2F("H_FD_neutral_beta_mom", "H_FD_neutral_beta_mom", 100, 0, 10.6, 100, 0, 1.2);
		H_FD_neutral_beta_mom.setTitle("NEUTRAL  #beta vs mom");
		H_FD_neutral_beta_mom.setTitleX("p (GeV)");
		H_FD_neutral_beta_mom.setTitleY("FTB #beta");
		H_FD_pos_mass_mom = new H2F("H_FD_pos_mass_mom", "H_FD_pos_mass_mom", 100, 0, 5, 150, -0.5, 2);
		H_FD_pos_mass_mom.setTitle("POS Mass^2 vs mom S%d");
		H_FD_pos_mass_mom.setTitleX("p (GeV)");
		H_FD_pos_mass_mom.setTitleY("M^2 (GeV^2)");
		H_FD_neg_mass_mom = new H2F("H_FD_neg_mass_mom", "H_FD_neg_mass_mom", 100, 0, 5, 150, -0.5, 2);
		H_FD_neg_mass_mom.setTitle("NEG Mass^2 vs mom S%d");
		H_FD_neg_mass_mom.setTitleX("p (GeV)");
		H_FD_neg_mass_mom.setTitleY("M^2 (GeV^2)");
		H_FD_neutral_mass_mom = new H2F("H_FD_neutral_mass_mom", "H_FD_neutral_mass_mom", 100, 0, 5, 150, -0.5, 2);
		H_FD_neutral_mass_mom.setTitle("NEUTRAL Mass^2 vs mom S%d");
		H_FD_neutral_mass_mom.setTitleX("p (GeV)");
		H_FD_neutral_mass_mom.setTitleY("M^2 (GeV^2)");
		H_FD_pos_mass_the = new H2F("H_FD_pos_mass_the", "H_FD_pos_mass_the", 100, 0, 45, 100, -0.5, 2);
		H_FD_pos_mass_the.setTitle("POS Mass^2 vs #theta");
		H_FD_pos_mass_the.setTitleX("#theta (^o)");
		H_FD_pos_mass_the.setTitleY("M^2 (GeV^2");
		H_FD_neg_mass_the = new H2F("H_FD_neg_mass_the", "H_FD_neg_mass_the", 100, 0, 45, 100, -0.5, 2);
		H_FD_neg_mass_the.setTitle("NEG Mass^2 vs #theta");
		H_FD_neg_mass_the.setTitleX("#theta (^o)");
		H_FD_neg_mass_the.setTitleY("M^2 (GeV^2");
		H_FD_neutral_mass_the = new H2F("H_FD_neutral_mass_the", "H_FD_neutral_mass_the", 100, 0, 45, 100, -0.5, 2);
		H_FD_neutral_mass_the.setTitle("NEUTRAL Mass^2 vs #theta");
		H_FD_neutral_mass_the.setTitleX("#theta (^o)");
		H_FD_neutral_mass_the.setTitleY("M^2 (GeV^2");
		H_FD_pos_mass_phi = new H2F("H_FD_pos_mass_phi", "H_FD_pos_mass_phi", 100, -180, 180, 100, -0.5, 2);
		H_FD_pos_mass_phi.setTitle("POS Mass^2 vs #phita");
		H_FD_pos_mass_phi.setTitleX("#phi (^o)");
		H_FD_pos_mass_phi.setTitleY("M^2 (GeV^2");
		H_FD_neg_mass_phi = new H2F("H_FD_neg_mass_phi", "H_FD_neg_mass_phi", 100, -180, 180, 100, -0.5, 2);
		H_FD_neg_mass_phi.setTitle("NEG Mass^2 vs #phita");
		H_FD_neg_mass_phi.setTitleX("#phi (^o)");
		H_FD_neg_mass_phi.setTitleY("M^2 (GeV^2");
		H_FD_neutral_mass_phi = new H2F("H_FD_neutral_mass_phi", "H_FD_neutral_mass_phi", 100, -180, 180, 100, -0.5, 2);
		H_FD_neutral_mass_phi.setTitle("NEUTRAL Mass^2 vs #phita");
		H_FD_neutral_mass_phi.setTitleX("#phi (^o)");
		H_FD_neutral_mass_phi.setTitleY("M^2 (GeV^2");
		
		// for calculated FTOF mass from paddle 1a, 1b for different sectors
		
		H_FTOF_pos_beta_mom_pad1a = new H2F[6];
		H_FTOF_neg_beta_mom_pad1a = new H2F[6];
		H_FTOF_pos_beta_mom_pad1b = new H2F[6];
		H_FTOF_neg_beta_mom_pad1b = new H2F[6];
		H_FTOF_pos_mass_mom_pad1a = new H2F[6];
		H_FTOF_pos_mass_the_pad1a = new H2F[6];
		H_FTOF_neg_mass_mom_pad1a = new H2F[6];
		H_FTOF_neg_mass_the_pad1a = new H2F[6];
		H_FTOF_pos_mass_mom_pad1b = new H2F[6];
		H_FTOF_pos_mass_the_pad1b = new H2F[6];
		H_FTOF_neg_mass_mom_pad1b = new H2F[6];
		H_FTOF_neg_mass_the_pad1b = new H2F[6];
		
		// mass from EB
		
		for (int s = 0; s < 6; s++) {
			
			H_FTOF_pos_beta_mom_pad1a[s] = new H2F(String.format("H_FTOF_pos_beta_mom_pad1a_%d", s + 1),
					String.format("H_FTOF_pos_beta_mom_pad1a_%d", s + 1), 100, 0, 10.6, 100, 0, 1.2);
			H_FTOF_pos_beta_mom_pad1a[s].setTitle(String.format("POS TOF1A #beta vs mom S%d", s + 1));
			H_FTOF_pos_beta_mom_pad1a[s].setTitleX("p (GeV)");
			H_FTOF_pos_beta_mom_pad1a[s].setTitleY("TOF #beta");
			H_FTOF_neg_beta_mom_pad1a[s] = new H2F(String.format("H_FTOF_pos_beta_neg_pad1a_%d", s + 1),
					String.format("H_FTOF_pos_beta_neg_pad1a_%d", s + 1), 100, 0, 10.6, 100, 0, 1.2);
			H_FTOF_neg_beta_mom_pad1a[s].setTitle(String.format("NEG TOF1A #beta vs mom S%d", s + 1));
			H_FTOF_neg_beta_mom_pad1a[s].setTitleX("p (GeV)");
			H_FTOF_neg_beta_mom_pad1a[s].setTitleY("TOF #beta");
			H_FTOF_pos_beta_mom_pad1b[s] = new H2F(String.format("H_FTOF_pos_beta_mom_pad1b_%d", s + 1),
					String.format("H_FTOF_pos_beta_mom_pad1b_%d", s + 1), 100, 0, 10.6, 100, 0, 1.2);
			H_FTOF_pos_beta_mom_pad1b[s].setTitle(String.format("POS TOF1B #beta vs mom S%d", s + 1));
			H_FTOF_pos_beta_mom_pad1b[s].setTitleX("p (GeV)");
			H_FTOF_pos_beta_mom_pad1b[s].setTitleY("TOF #beta");
			H_FTOF_neg_beta_mom_pad1b[s] = new H2F(String.format("H_FTOF_pos_beta_neg_pad1b_%d", s + 1),
					String.format("H_FTOF_pos_beta_neg_pad1b_%d", s + 1), 100, 0, 10.6, 100, 0, 1.2);
			H_FTOF_neg_beta_mom_pad1b[s].setTitle(String.format("NEG TOF1B #beta vs mom S%d", s + 1));
			H_FTOF_neg_beta_mom_pad1b[s].setTitleX("p (GeV)");
			H_FTOF_neg_beta_mom_pad1b[s].setTitleY("TOF #beta");
			H_FTOF_pos_mass_mom_pad1a[s] = new H2F(String.format("H_FTOF_pos_mass_mom_pad1a_%d", s + 1),
					String.format("H_FTOF_pos_mass_mom_pad1a_%d", s + 1), 100, 0, 5, 150, -0.5, 2);
			H_FTOF_pos_mass_mom_pad1a[s].setTitle(String.format("POS Mass^2 vs mom S%d", s + 1));
			H_FTOF_pos_mass_mom_pad1a[s].setTitleX("p (GeV)");
			H_FTOF_pos_mass_mom_pad1a[s].setTitleY("M^2 (GeV^2)");
			H_FTOF_pos_mass_the_pad1a[s] = new H2F(String.format("H_FTOF_pos_mass_the_pad1a_%d", s + 1),
					String.format("H_FTOF_pos_mass_the_pad1a_%d", s + 1), 100, 0, 45, 100, -0.5, 2);
			H_FTOF_pos_mass_the_pad1a[s].setTitle(String.format("POS Mass^2 vs #theta S%d", s + 1));
			H_FTOF_pos_mass_the_pad1a[s].setTitleX("#theta (^o)");
			H_FTOF_pos_mass_the_pad1a[s].setTitleY("M^2 (GeV^2)");
			H_FTOF_neg_mass_mom_pad1a[s] = new H2F(String.format("H_FTOF_neg_mass_mom_pad1a_%d", s + 1),
					String.format("H_FTOF_neg_mass_mom_pad1a_%d", s + 1), 100, 0, 5, 100, -0.5, 2);
			H_FTOF_neg_mass_mom_pad1a[s].setTitle(String.format("NEG Mass^2 vs mom S%d", s + 1));
			H_FTOF_neg_mass_mom_pad1a[s].setTitleX("p (GeV)");
			H_FTOF_neg_mass_mom_pad1a[s].setTitleY("M^2 (GeV^2)");
			H_FTOF_neg_mass_the_pad1a[s] = new H2F(String.format("H_FTOF_neg_mass_the_pad1a_%d", s + 1),
					String.format("H_FTOF_neg_mass_the_pad1a_%d", s + 1), 100, 0, 45, 100, -0.5, 2);
			H_FTOF_neg_mass_the_pad1a[s].setTitle(String.format("NEG Mass^2 vs #theta S%d", s + 1));
			H_FTOF_neg_mass_the_pad1a[s].setTitleX("#theta (^o)");
			H_FTOF_neg_mass_the_pad1a[s].setTitleY("M^2 (GeV^2)");
			H_FTOF_pos_mass_mom_pad1b[s] = new H2F(String.format("H_FTOF_pos_mass_mom_pad1b_%d", s + 1),
					String.format("H_FTOF_pos_mass_mom_pad1b_%d", s + 1), 100, 0, 5, 100, -0.5, 2);
			H_FTOF_pos_mass_mom_pad1b[s].setTitle(String.format("POS Mass^2 vs mom S%d", s + 1));
			H_FTOF_pos_mass_mom_pad1b[s].setTitleX("p (GeV)");
			H_FTOF_pos_mass_mom_pad1b[s].setTitleY("M^2 (GeV^2)");
			H_FTOF_pos_mass_the_pad1b[s] = new H2F(String.format("H_FTOF_pos_mass_the_pad1b_%d", s + 1),
					String.format("H_FTOF_pos_mass_the_pad1b_%d", s + 1), 100, 0, 45, 100, -0.5, 2);
			H_FTOF_pos_mass_the_pad1b[s].setTitle(String.format("POS Mass^2 vs #theta S%d", s + 1));
			H_FTOF_pos_mass_the_pad1b[s].setTitleX("#theta (^o)");
			H_FTOF_pos_mass_the_pad1b[s].setTitleY("M^2 (GeV^2)");
			H_FTOF_neg_mass_mom_pad1b[s] = new H2F(String.format("H_FTOF_neg_mass_mom_pad1b_%d", s + 1),
					String.format("H_FTOF_neg_mass_mom_pad1b_%d", s + 1), 100, 0, 5, 100, -0.5, 2);
			H_FTOF_neg_mass_mom_pad1b[s].setTitle(String.format("NEG Mass^2 vs mom S%d", s + 1));
			H_FTOF_neg_mass_mom_pad1b[s].setTitleX("p (GeV)");
			H_FTOF_neg_mass_mom_pad1b[s].setTitleY("M^2 (GeV^2)");
			H_FTOF_neg_mass_the_pad1b[s] = new H2F(String.format("H_FTOF_neg_mass_the_pad1b_%d", s + 1),
					String.format("H_FTOF_neg_mass_the_pad1b_%d", s + 1), 100, 0, 45, 100, -0.5, 2);
			H_FTOF_neg_mass_the_pad1b[s].setTitle(String.format("NEG Mass^2 vs #theta S%d", s + 1));
			H_FTOF_neg_mass_the_pad1b[s].setTitleX("#theta (^o)");
			H_FTOF_neg_mass_the_pad1b[s].setTitleY("M^2 (GeV^2)");
			
		}
		
		
	}
	
	
	public void fillRecBank(DataBank recBank) {
		STT = recBank.getFloat("startTime", 0);
		//RFT = recBank.getFloat("RFTime", 0);
	}	
	
//	
//	public void fillRecFTBank(DataBank recBank) {
//		
//		STT = recBank.getFloat("startTime", 0);
////		RFT = recBank.getFloat("RFTime", 0);		
//	}

	public void fillFTOF(DataBank part, DataBank bank) {
		for (int r = 0; r < bank.rows(); r++) {
			if (bank.getByte("detector", r) == 12) {
				
				if (bank.getShort("pindex", r) == pip_part_ind && bank.getByte("layer", r) == 2) {
					pip_FTOF_pad1b = bank.getShort("component", r);
					pip_FTOF1b_t = bank.getFloat("time", r);
					pip_FTOF1b_path = bank.getFloat("path", r);
					float pip_beta = pip_mom / (float) Math.sqrt(pip_mom * pip_mom + 0.13957f * 0.13957f);
					pip_FTOF1b_vt = pip_FTOF1b_t - pip_FTOF1b_path / (pip_beta * 29.98f) - STT;
					
				} // FTOF panal 1b 
				
				if (bank.getShort("pindex", r) > -1 && bank.getShort("pindex", r) < part.rows()) {
					byte q = part.getByte("charge", bank.getShort("pindex", r));
					float px = part.getFloat("px", bank.getShort("pindex", r));
					float py = part.getFloat("py", bank.getShort("pindex", r));
					float pz = part.getFloat("pz", bank.getShort("pindex", r));
					double mom = Math.sqrt(px * px + py * py + pz * pz);
					double the = Math.toDegrees(Math.acos(pz / mom));
					double TOFbeta = bank.getFloat("path", r) / (29.98f * (bank.getFloat("time", r) - STT));
					// double TOFmass = mom * Math.sqrt( 1/(TOFbeta*TOFbeta) - 1);
					double TOFmass = mom * mom * (1 / (TOFbeta * TOFbeta) - 1);
					int s = bank.getInt("sector", r) - 1;
					if (bank.getByte("layer", r) == 1 && q > 0) {
						H_FTOF_pos_beta_mom_pad1a[s].fill(mom, TOFbeta);
						H_FTOF_pos_mass_mom_pad1a[s].fill(mom, TOFmass);
						H_FTOF_pos_mass_the_pad1a[s].fill(the, TOFmass);
					}
					if (bank.getByte("layer", r) == 1 && q < 0) {
						H_FTOF_neg_beta_mom_pad1a[s].fill(mom, TOFbeta);
						H_FTOF_neg_mass_mom_pad1a[s].fill(mom, TOFmass);
						H_FTOF_neg_mass_the_pad1a[s].fill(the, TOFmass);
					}
					if (bank.getByte("layer", r) == 2 && q > 0) {
						H_FTOF_pos_beta_mom_pad1b[s].fill(mom, TOFbeta);
						H_FTOF_pos_mass_mom_pad1b[s].fill(mom, TOFmass);
						H_FTOF_pos_mass_the_pad1b[s].fill(the, TOFmass);
					}
					if (bank.getByte("layer", r) == 2 && q < 0) {
						H_FTOF_neg_beta_mom_pad1b[s].fill(mom, TOFbeta);
						H_FTOF_neg_mass_mom_pad1b[s].fill(mom, TOFmass);
						H_FTOF_neg_mass_the_pad1b[s].fill(the, TOFmass);
					}
				}

				
			} /// FTOF 
			
		}
	}
	
	
	public int makeFTElectron(DataBank bank) {
		int NFTElec = 0;
		//Mp = 0.93827f;
		for (int k = 0; k < bank.rows(); k++) {
			int pid = bank.getInt("pid", k);
			byte q = bank.getByte("charge", k);
			float px = bank.getFloat("px", k);
			float py = bank.getFloat("py", k);
			float pz = bank.getFloat("pz", k);
			int status = bank.getShort("status", k);
			boolean inFT = (status >= 1000 && status < 2000);
			e_mom = (float) Math.sqrt(px * px + py * py + pz * pz);
			if (pid == 11 && q == -1 && inFT && e_mom > 1.0 ) NFTElec++;
			
			if (NFTElec == 1) {
				e_the = (float) Math.toDegrees(Math.acos(pz / e_mom));
				e_phi = (float) Math.toDegrees(Math.atan2(py, px));
				//e_vx = bank.getFloat("vx", k);
				//e_vy = bank.getFloat("vy", k);
				Ve = new LorentzVector(px, py, pz, e_mom);
				VGS = new LorentzVector(0, 0, 0, 0);
				VGS.add(VB);
				VGS.sub(Ve);
				e_Q2 = (float) -VGS.mass2();
				e_xB = e_Q2 / (2f * Mp * (Eb - e_mom));
				e_W = (float) Math.sqrt(Mp * Mp + e_Q2 * (1f / e_xB - 1f));
				return k;
				
			}
		}
			
		
		return -1;
	}
	
	public void makeOthers(DataBank recbank,  DataBank recFTbank) {
		
		for (int k = 0; k < recbank.rows(); k++) {
			byte q = recbank.getByte("charge", k);
			int pid = recbank.getInt("pid", k);
			int ftbpid = recFTbank.getInt("pid", k);
			float px = recbank.getFloat("px", k);
			float py = recbank.getFloat("py", k);
			float pz = recbank.getFloat("pz", k);
			float vx = recbank.getFloat("vx", k);
			float vy = recbank.getFloat("vy", k);
			float vz = recbank.getFloat("vz", k);
			float be = recbank.getFloat("beta", k);
			float chi2pid = recFTbank.getFloat("chi2pid", k);
			float ftbbe = recFTbank.getFloat("beta", k);
			int status = recbank.getShort("status", k);
			boolean inDC = (status >= 2000 && status < 4000);
			float mom = (float) Math.sqrt(px * px + py * py + pz * pz);
			float the = (float) Math.toDegrees(Math.acos(pz / mom));
			float phi = (float) Math.toDegrees(Math.atan2(py, px));
			double TOFmass = mom * mom * (1 / ( ftbbe * ftbbe ) - 1);
			if (ftbpid == 211 && pip_part_ind == -1 && inDC && mom > 1.5 ) {
				pip_part_ind = k;
				pip_mom = mom;
				pip_the = the;
				pip_phi = (float) Math.toDegrees(Math.atan2(py, px));
				pip_vx = vx;
				pip_vy = vy;
				pip_vz = vz;
				pip_ftb_beta = ftbbe;
				Vpip = new LorentzVector(px, py, pz, Math.sqrt(pip_mom * pip_mom + 0.13957f * 0.13957f));

			}
			
			if (ftbpid == 2212 && prot_part_ind == -1 && inDC && mom > 1.5) {
				prot_part_ind = k;
				prot_mom = mom;
				prot_the = the;
				prot_phi = (float) Math.toDegrees(Math.atan2(py, px));
				prot_vx = vx;
				prot_vy = vy;
				prot_vz = vz;
				prot_ftb_beta = ftbbe;
				Vprot = new LorentzVector(px, py, pz, Math.sqrt(prot_mom * prot_mom + 0.93827f * 0.93827f));
				
			}
			
			if ( q > 0 && inDC && Math.abs(chi2pid) < 5.0) {
				H_FD_pos_beta_mom.fill(mom, ftbbe);
				H_FD_pos_mass_mom.fill(mom, TOFmass);
				H_FD_pos_mass_the.fill(the, TOFmass);
				H_FD_pos_mass_phi.fill(phi, TOFmass);
			}
			if ( q < 0 && inDC && Math.abs(chi2pid) < 5.0) {
				H_FD_neg_beta_mom.fill(mom, ftbbe);
				H_FD_neg_mass_mom.fill(mom, TOFmass);
				H_FD_neg_mass_the.fill(the, TOFmass);
				H_FD_neg_mass_phi.fill(phi, TOFmass);

			}
			if (q == 0 && inDC && Math.abs(chi2pid) < 5.0) {
				H_FD_neutral_beta_mom.fill(mom, ftbbe);
				H_FD_neutral_mass_mom.fill(mom, TOFmass);
				H_FD_neutral_mass_the.fill(the, TOFmass);
				H_FD_neutral_mass_phi.fill(phi, TOFmass);
			}
			
		} // FOR LOOP
	
	} //MAKEOTHER
	
	public void resetCounters() {
		e_ft_part_ind = -1;
		pip_part_ind = -1;
		prot_part_ind = -1;
		pip_FTOF_pad1b = -1;
			
	}
	
	

		
	
	public void processEvent(DataEvent event) {
		resetCounters();
		if (event.hasBank("RECFT::Event"))
			fillRecBank(event.getBank("RECFT::Event"));
//		if (event.hasBank("REC::Event"))
//			fillRecBank(event.getBank("REC::Event"));
		if (event.hasBank("REC::Particle"))
			e_ft_part_ind = makeFTElectron(event.getBank("REC::Particle"));
		if(e_ft_part_ind > -1) {
			
			if (event.hasBank("REC::Particle") && event.hasBank("RECFT::Particle")) makeOthers(event.getBank("REC::Particle"), event.getBank("RECFT::Particle"));
			if (event.hasBank("REC::Scintillator")) fillFTOF(event.getBank("REC::Particle"), event.getBank("REC::Scintillator"));
			
			
			FillHists();
			
		} // e_ft_part_ind > -1
		
		
	} //processEvent
	
	
	
	public void FillHists() {
		
		H_FT_e_t_f.fill(e_phi, e_the);
		H_FT_e_p_f.fill(e_phi, e_mom);
		H_FT_e_p_the.fill(e_the, e_mom);
		H_FT_W_Q2.fill(e_W, e_Q2);
		H_FT_e_xB_Q2.fill(e_xB, e_Q2);
		H_FT_W.fill(e_W);
		
	}
	
	
	public void plot() {
		EmbeddedCanvas can_e_overview = new EmbeddedCanvas();
		can_e_overview.setSize(3600, 1800);
		can_e_overview.divide(3, 2);
		can_e_overview.setAxisTitleSize(24);
		can_e_overview.setAxisFontSize(24);
		can_e_overview.setTitleSize(24);
		can_e_overview.cd(0);
		can_e_overview.draw(H_FT_e_t_f);
		can_e_overview.cd(1);
		can_e_overview.draw(H_FT_e_p_f);
		can_e_overview.cd(2);
		can_e_overview.draw(H_FT_e_p_the);
		can_e_overview.cd(3);
		can_e_overview.draw(H_FT_W_Q2);
		can_e_overview.cd(4);
		can_e_overview.draw(H_FT_e_xB_Q2);
		can_e_overview.cd(5);
		can_e_overview.draw(H_FT_W);
		can_e_overview.getPad(5).getAxisZ().setRange(5, 25);
		can_e_overview.save(String.format("/home/akhanal/eclipse-workspace/myAnalyser/plots/dst_FT_e_overview.png"));
		System.out.println("saved plots dst_FT_e_overview.png");
		
		EmbeddedCanvas can_mass_overview = new EmbeddedCanvas();
		can_mass_overview.setSize(3600, 1800);
		can_mass_overview.divide(3, 2);
		can_mass_overview.setAxisTitleSize(24);
		can_mass_overview.setAxisFontSize(24);
		can_mass_overview.setTitleSize(24);
		can_mass_overview.cd(0);
		can_mass_overview.draw(H_FD_pos_mass_mom);
		can_mass_overview.cd(1);
		can_mass_overview.draw(H_FD_neg_mass_mom);
		can_mass_overview.cd(2);
		can_mass_overview.draw(H_FD_neutral_mass_mom);
		can_mass_overview.cd(3);
		can_mass_overview.draw(H_FD_pos_mass_the);
		can_mass_overview.cd(4);
		can_mass_overview.draw(H_FD_neg_mass_the);
		can_mass_overview.cd(5);
		can_mass_overview.draw(H_FD_neutral_mass_the);
		can_mass_overview.save(String.format("/home/akhanal/eclipse-workspace/myAnalyser/plots/dst_mass_overview.png"));
		System.out.println("saved plots dst_mass_overview.png");
		
		EmbeddedCanvas can_massvsphi_overview = new EmbeddedCanvas();
		can_massvsphi_overview.setSize(1800, 1200);
		can_massvsphi_overview.divide(3, 1);
		can_massvsphi_overview.setAxisTitleSize(24);
		can_massvsphi_overview.setAxisFontSize(24);
		can_massvsphi_overview.setTitleSize(24);
		can_massvsphi_overview.cd(0);
		can_massvsphi_overview.draw(H_FD_pos_mass_phi);
		can_massvsphi_overview.cd(1);
		can_massvsphi_overview.draw(H_FD_neg_mass_phi);
		can_massvsphi_overview.cd(2);
		can_massvsphi_overview.draw(H_FD_neutral_mass_phi);
		can_massvsphi_overview.save(String.format("/home/akhanal/eclipse-workspace/myAnalyser/plots/dst_massvsphi_overview.png"));
		System.out.println("saved plots dst_massvsphi_overview.png");
		
		EmbeddedCanvas can_bevsp_overview = new EmbeddedCanvas();
		can_bevsp_overview.setSize(3600, 1800);
		can_bevsp_overview.divide(3, 2);
		can_bevsp_overview.setAxisTitleSize(24);
		can_bevsp_overview.setAxisFontSize(24);
		can_bevsp_overview.setTitleSize(24);
		can_bevsp_overview.cd(0);
		can_bevsp_overview.draw(H_FD_pos_beta_mom);	
		can_bevsp_overview.draw(F_prot_beta_mom, "same");
		can_bevsp_overview.draw(F_kp_beta_mom, "same");
		can_bevsp_overview.draw(F_pip_beta_mom, "same");
		can_bevsp_overview.cd(1);
		can_bevsp_overview.draw(H_FD_neg_beta_mom);
		can_bevsp_overview.draw(F_prot_beta_mom, "same");
		can_bevsp_overview.draw(F_kp_beta_mom, "same");
		can_bevsp_overview.draw(F_pip_beta_mom, "same");
		can_bevsp_overview.cd(2);
		can_bevsp_overview.draw(H_FD_neutral_beta_mom);
		H1F H_FD_pos_mass_mom_projY = H_FD_pos_mass_mom.projectionY();
		H_FD_pos_mass_mom_projY.setTitle("POS M^2 (GeV^2)");
		can_bevsp_overview.cd(3);
		can_bevsp_overview.draw(H_FD_pos_mass_mom_projY);
		H1F H_FD_neg_mass_mom_projY = H_FD_neg_mass_mom.projectionY();
		H_FD_neg_mass_mom_projY.setTitle("NEG M^2 (GeV^2)");
		can_bevsp_overview.cd(4);
		can_bevsp_overview.draw(H_FD_neg_mass_mom_projY);
		H1F H_FD_neutral_mass_mom_projY = H_FD_neutral_mass_mom.projectionY();
		H_FD_neutral_mass_mom_projY.setTitle("NEUTRAL M^2 (GeV^2)");
		can_bevsp_overview.cd(5);
		can_bevsp_overview.draw(H_FD_neutral_mass_mom_projY);
		//can_bevsp_overview.save(String.format("/home/akhanal/eclipse-workspace/rg-b/src/plots/dst_bevs_overview.png"));
		can_bevsp_overview.save(String.format("/home/akhanal/eclipse-workspace/myAnalyser/plots/dst_bevsp_overview.png"));
		System.out.println("saved plots dst_bevsp_overview.png");
		
		
		
		EmbeddedCanvas can_e_FTOF1A_mass = new EmbeddedCanvas();
		can_e_FTOF1A_mass.setSize(3600, 4800);
		can_e_FTOF1A_mass.divide(6, 8);
		can_e_FTOF1A_mass.setAxisTitleSize(24);
		can_e_FTOF1A_mass.setAxisFontSize(24);
		can_e_FTOF1A_mass.setTitleSize(24);
		for (int s = 0; s < 6; s++) {
			
			H1F H_FTOF_pos_mass_mom_pad1a_projY = H_FTOF_pos_mass_mom_pad1a[s].projectionY();
			H_FTOF_pos_mass_mom_pad1a_projY.setTitle(String.format("POS TOF1A mass S%d", s + 1));
			H_FTOF_pos_mass_mom_pad1a_projY.setTitleX("M^2 (GeV^2)");
			H1F H_FTOF_neg_mass_mom_pad1a_projY = H_FTOF_neg_mass_mom_pad1a[s].projectionY();
			H_FTOF_neg_mass_mom_pad1a_projY.setTitle(String.format("NEG TOF1A mass S%d", s + 1));
			H_FTOF_neg_mass_mom_pad1a_projY.setTitleX("M^2 (GeV^2)");
			can_e_FTOF1A_mass.cd(s);
			can_e_FTOF1A_mass.getPad(s).getAxisY().setLog(true);
			can_e_FTOF1A_mass.draw(H_FTOF_pos_mass_mom_pad1a_projY);
			can_e_FTOF1A_mass.cd(s + 6);
			can_e_FTOF1A_mass.draw(H_FTOF_pos_mass_mom_pad1a[s]);
			can_e_FTOF1A_mass.cd(s + 12);
			can_e_FTOF1A_mass.draw(H_FTOF_pos_mass_the_pad1a[s]);
			can_e_FTOF1A_mass.cd(s + 18);
			can_e_FTOF1A_mass.draw(H_FTOF_pos_beta_mom_pad1a[s]);
			can_e_FTOF1A_mass.cd(s + 24);
			can_e_FTOF1A_mass.draw(H_FTOF_neg_mass_mom_pad1a_projY);
			can_e_FTOF1A_mass.cd(s + 30);
			can_e_FTOF1A_mass.draw(H_FTOF_neg_mass_mom_pad1a[s]);
			can_e_FTOF1A_mass.cd(s + 36);
			can_e_FTOF1A_mass.draw(H_FTOF_neg_mass_the_pad1a[s]);
			can_e_FTOF1A_mass.cd(s + 42);
			can_e_FTOF1A_mass.draw(H_FTOF_neg_beta_mom_pad1a[s]);
			
			}
		
			can_e_FTOF1A_mass.save(String.format("/home/akhanal/eclipse-workspace/myAnalyser/plots/dst_FTOF1A_mass.png"));
			System.out.println(String.format("saved plots/dst_FTOF1A_mass.png"));
			
			EmbeddedCanvas can_e_FTOF1B_mass = new EmbeddedCanvas();
			can_e_FTOF1B_mass.setSize(3600, 4800);
			can_e_FTOF1B_mass.divide(6, 8);
			can_e_FTOF1B_mass.setAxisTitleSize(24);
			can_e_FTOF1B_mass.setAxisFontSize(24);
			can_e_FTOF1B_mass.setTitleSize(24);
			for (int s = 0; s < 6; s++) {
				H1F H_FTOF_pos_mass_mom_pad1b_projY = H_FTOF_pos_mass_mom_pad1b[s].projectionY();
				H_FTOF_pos_mass_mom_pad1b_projY.setTitle(String.format("POS TOF1B mass^2 S%d", s + 1));
				H_FTOF_pos_mass_mom_pad1b_projY.setTitleX("M^2 (GeV^2)");
				H1F H_FTOF_neg_mass_mom_pad1b_projY = H_FTOF_neg_mass_mom_pad1b[s].projectionY();
				H_FTOF_neg_mass_mom_pad1b_projY.setTitle(String.format("NEG TOF1B mass^2 S%d", s + 1));
				H_FTOF_neg_mass_mom_pad1b_projY.setTitleX("M^2 (GeV^2)");
				can_e_FTOF1B_mass.cd(s);
				can_e_FTOF1B_mass.getPad(s).getAxisY().setLog(true);
				can_e_FTOF1B_mass.draw(H_FTOF_pos_mass_mom_pad1b_projY);
				can_e_FTOF1B_mass.cd(s + 6);
				can_e_FTOF1B_mass.draw(H_FTOF_pos_mass_mom_pad1b[s]);
				can_e_FTOF1B_mass.cd(s + 12);
				can_e_FTOF1B_mass.draw(H_FTOF_pos_mass_the_pad1b[s]);
				can_e_FTOF1B_mass.cd(s + 18);
				can_e_FTOF1B_mass.draw(H_FTOF_pos_beta_mom_pad1b[s]);
				can_e_FTOF1B_mass.cd(s + 24);
				can_e_FTOF1B_mass.draw(H_FTOF_neg_mass_mom_pad1b_projY);
				can_e_FTOF1B_mass.cd(s + 30);
				can_e_FTOF1B_mass.draw(H_FTOF_neg_mass_mom_pad1b[s]);
				can_e_FTOF1B_mass.cd(s + 36);
				can_e_FTOF1B_mass.draw(H_FTOF_neg_mass_the_pad1b[s]);
				can_e_FTOF1B_mass.cd(s + 42);
				can_e_FTOF1B_mass.draw(H_FTOF_neg_beta_mom_pad1b[s]);
			}
			
			can_e_FTOF1B_mass.save(String.format("/home/akhanal/eclipse-workspace/myAnalyser/plots/dst_FTOF1B_mass.png"));
			System.out.println(String.format("saved plots/dst_FTOF1B_mass.png"));
				
	}
	
	public void write() {
		
		TDirectory dirout = new TDirectory();
		dirout.mkdir("/FTElec/");
		dirout.cd("/FTElec/");
		dirout.addDataSet(H_FT_e_t_f, H_FT_e_p_f, H_FT_W_Q2);
	}
	
	

	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	public static void main(String[] args) {
		System.setProperty("java.awt.headless", "true");
		GStyle.setPalette("kRainBow");
		int count = 0;
		int maxevents = 5000000;
		ft_ana ana = new ft_ana();
		File fileIn = new File(args[0]);
		java.util.Date date1 = new java.util.Date();
		System.out.println(date1);
		HipoDataSource reader = new HipoDataSource();
		reader.open(fileIn);
		while (reader.hasEvent() && count < maxevents) {
			DataEvent event = reader.getNextEvent();
			ana.processEvent(event);
			count++;
						
			
		}
		
		System.out.println("Total events : " + count);
		reader.close();
		ana.plot();
		ana.write();
		System.out.println("Bye.");
	}
		
	
	
}
