package ImageAnalysis;

import java.io.File;
import java.util.List;

public class Runner {

	public static void main(String[] args) {



		DataPreparation dataPreparation = new DataPreparation();
		ImageCreator imageCreator = new ImageCreator();
		LBP lbp = new LBP();
		HistogramCreator hist = new HistogramCreator();

		List<UserWalk> usersWalks = dataPreparation.getUsersWalks();
		dataPreparation.generateMinMax(usersWalks);

		for(UserWalk userWalk : usersWalks) {

			String subj = "subj_";

			if (userWalk.getUser() / 100 >= 1) {
				subj += String.valueOf(userWalk.getUser());
			} else if (userWalk.getUser() / 10 >= 1) {
				subj += "0" + userWalk.getUser();
			} else {
				subj += "00" + userWalk.getUser();
			}

			String path = subj + File.pathSeparator+"walk" + userWalk.getWalk();

			// crea i file con i dati interpolati
			dataPreparation.dataInterpolation(Constants.datasetSource+File.pathSeparator+path, userWalk);

			// crea i file con i dati per la creazione di ogni singolo pixel
			dataPreparation.createImageDataset(Constants.datasetSource+File.pathSeparator+path);

			// crea l'immagine
			imageCreator.createImage(Math.ceil(dataPreparation.getMeanStepLength()), Constants.datasetSource + subj +File.pathSeparator+"walk" + userWalk.getWalk(), Constants.imagesSource + subj+File.pathSeparator+ userWalk.getWalk()+File.pathSeparator, userWalk);

			// LBP
			lbp.calculateLBP(Constants.imagesSource + subj +File.pathSeparator+ userWalk.getWalk() +File.pathSeparator+ userWalk.getWalk(), userWalk, Math.ceil(dataPreparation.getMeanStepLength()));

			// crea gli istogrammi
			hist.createHistogramUserWalk(userWalk, Constants.imagesSource + subj+File.pathSeparator+userWalk.getWalk()+File.pathSeparator);
		}

		// compara con chi-square e bhattacharyya
		HistogramComparator histogramComparator = new HistogramComparator();
		histogramComparator.compare(usersWalks);

		// compara con correlation e intersection
		HistogramComparatorInverted histogramComparatorInverted = new HistogramComparatorInverted();
		histogramComparatorInverted.compare(usersWalks);
	}
}