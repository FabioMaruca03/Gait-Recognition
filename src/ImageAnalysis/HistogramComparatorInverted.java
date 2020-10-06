package ImageAnalysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.stream.Collectors;

import org.opencv.imgproc.Imgproc;

/**
 * questo comparatore viene utilizzato solo per correlation e intersection
 * siccome viene calcolata la somiglianza e non la distanza tra istogrammi, le funzioni sono le stesse presenti in HistogramComparator
 * ma le funzioni vengono adattate al comportamento delle metriche
 * 
 */
public class HistogramComparatorInverted {
	
	List<UserWalk> gallery = new ArrayList<>();
	LinkedHashMap<UserWalk, List<Double>> probe = new LinkedHashMap<>();
	LinkedHashMap<UserWalk, List<Double>> probeCollected = new LinkedHashMap<>();
	Double min = Double.MAX_VALUE;
	Double max = 0.0;
	
	long totalGA = 0;
	long totalFA = 0;
	long totalGR = 0;
	long totalFR = 0;
	
	int localGA = 0;
	int localFA = 0;
	int localGR = 0;
	int localFR = 0;
	
	List<Double> far = new ArrayList<>();
	List<Double> frr = new ArrayList<>();
	List<Double> roc = new ArrayList<>();

	public void compare(List<UserWalk> usersWalks) {
		for(int i = 0; i < 50; i++) {
			for(UserWalk userWalk : usersWalks) {
				if(userWalk.getUser() == i) {
					if (userWalk.getWalk() < 6) {
						gallery.add(userWalk);
					} else {
						probe.put(userWalk, new ArrayList<>());
					}
				}
			}
		}
		
		for (UserWalk probeWalk : probe.keySet()) {
			// inizializza la lista delle distanze rispetto alle immagini della gallery
			List<Double> distance = new ArrayList<>();
			
			for (UserWalk galleryWalk : gallery) {
				// aggiunge la distanza alla lista
				distance.add(calculateHistDistance(probeWalk, galleryWalk));
			}

			probe.put(probeWalk, distance);
		}

		calculateCmc();
		findMinMax();
		normalizeResultProbe();
		collectResults();
		normalizeResultProbeCollected();
		verify();
		getClosestResult();
	}

	private void findMinMax() {
		for (UserWalk userWalk : probe.keySet()) {

			for (int i = 0; i < probe.get(userWalk).size(); i++) {
				// setta il massimo
				max = Math.max(max, probe.get(userWalk).get(i));
				// setta il minimo
				min = Math.min(min, probe.get(userWalk).get(i));
			}
		}
	}
	
	private Double calculateHistDistance(UserWalk probeWalk, UserWalk galleryWalk) {

		// SOLO 1 E 3
//		diffR = Imgproc.compareHist(probeWalk.getHistograms().get(0), galleryWalk.getHistograms().get(0), 1);
//		diffG = Imgproc.compareHist(probeWalk.getHistograms().get(1), galleryWalk.getHistograms().get(1), 1);
//		diffB = Imgproc.compareHist(probeWalk.getHistograms().get(2), galleryWalk.getHistograms().get(2), 1);
		
		double diffR = Imgproc.compareHist(probeWalk.getHistograms().get(0), galleryWalk.getHistograms().get(0), 3);
		double diffG = Imgproc.compareHist(probeWalk.getHistograms().get(1), galleryWalk.getHistograms().get(1), 3);
		double diffB = Imgproc.compareHist(probeWalk.getHistograms().get(2), galleryWalk.getHistograms().get(2), 3);
		
//		return (diffR + diffG + diffB)/3;
		return ((0.2*diffR) + (0.7*diffG) + (0.1*diffB));
	}

	private void collectResults() {
		for (UserWalk userWalk : probe.keySet()) {
			List<Double> collected = new ArrayList<>();
			List<Double> tempDepo = new ArrayList<>();

			for (int i = 0; i < probe.get(userWalk).size() + 1; i++) {
				if (i % 6 != 0 || i == 0) {
					tempDepo.add(probe.get(userWalk).get(i));
				} else {
					Collections.sort(tempDepo);
					collected.add(tempDepo.get(tempDepo.size()-1));
					
					max = Math.max(max, tempDepo.get(tempDepo.size()-1));

					min = Math.min(min, tempDepo.get(tempDepo.size()-1));

					tempDepo.clear();
					if (i != probe.get(userWalk).size()) {
						tempDepo.add(probe.get(userWalk).get(i));
					}
				}
			}
			probeCollected.put(userWalk, collected);
		}
	}	

	private void normalizeResultProbe() {
		
		for (UserWalk userWalk : probe.keySet()) {
			List<Double> normalizedValues = new ArrayList<>();
			for (Double d : probe.get(userWalk)) {
				normalizedValues.add((d - min) / (max - min));
			}
			probe.replace(userWalk, normalizedValues);
		}
	}

	private void normalizeResultProbeCollected() {

		for (UserWalk userWalk : probeCollected.keySet()) {
			List<Double> normalizedValues = new ArrayList<>();
			for (Double d : probeCollected.get(userWalk)) {
				normalizedValues.add((d - min) / (max - min));
			}
			probeCollected.replace(userWalk, normalizedValues);
		}
	}
	
	private void getClosestResult() {

		int totalSuccess = 0;
		int totalFail = 0;

		for (UserWalk userWalk : probe.keySet()) {
			
			int indexOfMinimum = probe.get(userWalk).indexOf(Collections.max(probe.get(userWalk)));
			
			if (userWalk.getUser() == gallery.get(indexOfMinimum).getUser()) {
				totalSuccess++;
			} else {
				totalFail++;
			}
		}
	}
	

	public void calculateCmc() {

		List<Double> cmc = new ArrayList<>();

		for (int rank = 1; rank < gallery.size() + 1; rank++) {
			double cmsk = 0.0;
			double success = 0.0;

			for (UserWalk userWalk : probe.keySet()) {
				List<DiffObject> diffObjList = new ArrayList<>();

				// crea una copia del probe in cui viene salvato il valore con il relativo indice
				for (int i = 0; i < probe.get(userWalk).size(); i++) {
					diffObjList.add(new DiffObject(probe.get(userWalk).get(i), i));
				}

				diffObjList.sort(Collections.reverseOrder(Comparator.comparing(DiffObject::getDiff)));

				if (diffObjList.stream().limit(rank)
						.anyMatch(e -> gallery.get(e.index).getUser() == userWalk.getUser())) {
					success++;
				}
			}

			cmsk = success/probe.keySet().size();

			cmc.add(cmsk);
		}
	}

	private void verify() {
		
		for (int t = 10000; t > 0 ; t--) {
			Double threshold = t / 10000.0;
			localGA = 0;
			localFA = 0;
			localGR = 0;
			localFR = 0;
			
			for (UserWalk userWalk : probeCollected.keySet()) {
				for (int i = 0; i < probeCollected.get(userWalk).size(); i++) {
					checkIdentityClaim(userWalk, i, threshold);
				}
			}
			
			far.add(localFA/14112.0);
			frr.add(localFR/294.0);
			
			totalGA += localGA;
			totalFA += localFA;
			totalGR += localGR;
			totalFR += localFR;
		}
		
		roc = frr.stream().map(e -> 1 - e).collect(Collectors.toList());
	}
	
	private void checkIdentityClaim(UserWalk userWalk, int identity, double threshold) {
		
		if(probeCollected.get(userWalk).get(identity) >= threshold) {
			if(userWalk.getUser() == identity + 1) {
				localGA++;
			} else {
				localFA++;
			}
		} else {
			if(userWalk.getUser() == identity + 1) {
				localFR++;
			} else {
				localGR++;
			}
		}
	}
	
	
	class DiffObject {
		Double diff;
		int index;

		public DiffObject(Double diff, int index) {
			this.diff = diff;
			this.index = index;
		}

		public Double getDiff() { return diff; }

		public void setDiff(Double diff) { this.diff = diff; }

		public int getIndex() { return index; }

		public void setIndex(int index) { this.index = index; }

		@Override
		public String toString() {
			String res = "diff: " + diff + ", index: " + index;
			return res;
		}
	}
}
