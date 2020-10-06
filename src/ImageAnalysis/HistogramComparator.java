package ImageAnalysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.stream.Collectors;

import org.opencv.imgproc.Imgproc;

/**
 * questo comparatore viene utilizzato solo per chi-square e bhattacharyya
 *
 */
public class HistogramComparator {

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

	// prepara le mappe di probe e gallery
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
		
		// confronta ogni probe con tutti i template della gallery
		for (UserWalk probeWalk : probe.keySet()) {
			// inizializza la lista delle distanze rispetto alle immagini della gallery
			List<Double> distance = new ArrayList<>();
			
			for (UserWalk galleryWalk : gallery) {
				// aggiunge la distanza alla lista
				distance.add(calculateHistDistance(probeWalk, galleryWalk));
			}

			// associa il risultato ad ogni probe
			probe.put(probeWalk, distance);
		}

		calculateCmc();
		// se si normalizza il probeCollected e' inutile, quindi si puo' commentare
		findMinMax();
		// se si normalizza il probeCollected e' inutile, quindi si puo' commentare
		normalizeResultProbe();
		collectResults();
		// se si normalizza il probe e' inutile, quindi si puo' commentare
		normalizeResultProbeCollected();
		verify();
		getClosestResult();
	}

	/**
	 * calcola minimo e massimo per la normalizzazione di tutto il probe intero
	 */
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
	
	/**
	 * calcola la distanza tra gli istogrammi delle due camminate selezionate
	 * 
	 * @param probeWalk è la camminata dell'utente preso dalla lista dei probe
	 * @param galleryWalk è la camminata dell'autente preso dalla gallery
	 * @return la media pesata delle distanze degli istogrammi
	 */
	private Double calculateHistDistance(UserWalk probeWalk, UserWalk galleryWalk) {
		
		// la distanza viene calcolata per tutti gli assi separatamente
		// l'istogramma 0 corrisponde al canale rosso, 1 al verde e 2 al blu 
		// la funzione compareHist di OpenCV prende in input gli istogrammi delle camminate e la metrica di confronto
		// la metrica 1 è Correlation, 2 Chi-Square, 3 Intersection e 4 Bhattacharyya
		// qua vengono applicati solamente chi-square e bhattacharyya
		
//		diffR = Imgproc.compareHist(probeWalk.getHistograms().get(0), galleryWalk.getHistograms().get(0), 2);
//		diffG = Imgproc.compareHist(probeWalk.getHistograms().get(1), galleryWalk.getHistograms().get(1), 2);
//		diffB = Imgproc.compareHist(probeWalk.getHistograms().get(2), galleryWalk.getHistograms().get(2), 2);

		double diffR = Imgproc.compareHist(probeWalk.getHistograms().get(0), galleryWalk.getHistograms().get(0), 4);
		double diffG = Imgproc.compareHist(probeWalk.getHistograms().get(1), galleryWalk.getHistograms().get(1), 4);
		double diffB = Imgproc.compareHist(probeWalk.getHistograms().get(2), galleryWalk.getHistograms().get(2), 4);
		
		// viene restituito il valore della media che puo' essere pesata o aritmetica
//		return (diffR + diffG + diffB)/3;
		return ((0.2*diffR) + (0.7*diffG) + (0.1*diffB));
	}

	/**
	 * prende i migliori risultati per ogni user della gallery
	 */
	private void collectResults() {
		// itera su tutti i probe
		for (UserWalk userWalk : probe.keySet()) {
			// inizializza due liste d'appoggio
			List<Double> collected = new ArrayList<>();
			List<Double> tempDepo = new ArrayList<>();

			// itera sulle distanze di ogni probe
			for (int i = 0; i < probe.get(userWalk).size() + 1; i++) {
				// in questo caso si sa gia' che ogni subject ha 6 template, in caso il numero sia differente questo valore va regolato di conseguenza
				if (i % 6 != 0 || i == 0) {
					tempDepo.add(probe.get(userWalk).get(i));
				} else {
					// ogni sesto elemento si ordina la lista riempita finora e aggiunge il minimo ad un'altra lista
					Collections.sort(tempDepo);
					collected.add(tempDepo.get(0));
					
					// massimo e minimo vengono usati solo se si sceglie di normalizzare il probeCollected
					// setta il massimo
					max = Math.max(max, tempDepo.get(0));

					// setta il minimo
					min = Math.min(min, tempDepo.get(0));

					// svuota la lista con i set di valori
					tempDepo.clear();
					
					// aggiunge il sesto elemento analizzato ad ogni ciclo che altrimenti resterebbe escluso
					if (i != probe.get(userWalk).size()) {
						tempDepo.add(probe.get(userWalk).get(i));
					}
				}
			}
			// associa la lista con i risultati migliori al probe
			probeCollected.put(userWalk, collected);
		}
	}	

	/**
	 * normalizza tutte le distanze dei probe
	 */
	private void normalizeResultProbe() {
		
		// itera su tutti i probe probe
		for (UserWalk userWalk : probe.keySet()) {
			// crea una lista di appoggio
			List<Double> normalizedValues = new ArrayList<>();
			for (Double d : probe.get(userWalk)) {
				// normalizza le distanze e le salva nella lista d'appoggio
				normalizedValues.add((d - min) / (max - min));
			}
			// sostituisce i risultati del probe
			probe.replace(userWalk, normalizedValues);
		}
	}

	/**
	 * normalizza le distanze del probe solo per i migliori risultati con ogni subject della gallery
	 */
	private void normalizeResultProbeCollected() {

		// itera sui probe minimizzati
		for (UserWalk userWalk : probeCollected.keySet()) {
			// crea una lista di appoggio
			List<Double> normalizedValues = new ArrayList<>();
			for (Double d : probeCollected.get(userWalk)) {
				// normalizza le distanze e le salva nella lista d'appoggio
				normalizedValues.add((d - min) / (max - min));
			}
			// sostituisce i risultati nel probeColleted
			probeCollected.replace(userWalk, normalizedValues);
		}
	}
	
	/**
	 * identificazione closed set
	 */
	private void getClosestResult() {

		// inizializza success e fail
		int totalSuccess = 0;
		int totalFail = 0;

		// itera su tutte camminate del probe
		for (UserWalk userWalk : probe.keySet()) {
			
			// prende l'indice della camminata con distanza minore
			int indexOfMinimum = probe.get(userWalk).indexOf(Collections.min(probe.get(userWalk)));
			
			// se il subject in corrispondenza dell'indice e' uguale al probe incrementa il success
			if (userWalk.getUser() == gallery.get(indexOfMinimum).getUser()) {
				totalSuccess++;
			} else {
				totalFail++;
			}
		}
	}
	
	/**
	 * calcola i valori di CMS per ogni rango a partire da 1 fino alla misura della lunghezza della gallery
	 * @return la lista di valori che compongono la curva di CMC 
	 */
	public List<Double> calculateCmc() {
		
		List<Double> cmc = new ArrayList<>();

		// itera sui ranghi
		for(int rank = 1; rank < gallery.size() + 1; rank++) {
			double successRate = 0.0;
			
			// itera sulla lista dei probe
			for(UserWalk userWalk : probe.keySet()) {
				
				// l'oggetto DiffObject contiene il valore della distanza calcolata e la posizione all'interno della lista delle distanze per ogni camminata
				// in questo modo è possibile lavorare su una copia del probe senza perdere l'informazione relativa alla corrispondenza con il subject della gallery
				List<DiffObject> diffObjList = new ArrayList<>();

				// crea una copia del probe in cui viene salvato il valore con il relativo indice
				for (int i = 0; i < probe.get(userWalk).size(); i++) {
					diffObjList.add(new DiffObject(probe.get(userWalk).get(i), i));
				}
				
				// ordina la lista per valori crescenti della distanza
				diffObjList.sort(Comparator.comparing(DiffObject::getDiff));
				
				//  se i primi k elementi della lista ordinata contengono almeno un match tra subject di probe e gallery incrementa il success rate
				if(diffObjList.stream().limit(rank).anyMatch(e -> gallery.get(e.index).getUser() == userWalk.getUser())) {
					successRate++;
				}
			}
			
			// aggiunge alla lista di CMC il CMS calcolato per il rango k
			cmc.add(successRate/probe.keySet().size());
		}

		// restituisce la lista dei CMS calcolati per tutti i ranghi
		return cmc;
	}
	
	// VERIFICA
	/**
	 * esegue il test di verifica sull'intero probe
	 */
	private void verify() {

		// la verifica viene eseguita per valori diversi del threshold
		for (int t = 0; t < 10000; t++) {
			// viene fissato il threshold con variazioni di un decimo di millesimo
			Double threshold = t / 10000.0;
			localGA = 0;
			localFA = 0;
			localGR = 0;
			localFR = 0;
			
			// vengono valutate le dichiarazioni di ogni probe
			// probeCollected contiene le distanze migliori relative ad ogni subject della gallery
			for (UserWalk userWalk : probeCollected.keySet()) {
				
				// l'identitià dichiarata dal probe cambia ad ogni iterazione
				for (int i = 1; i <= probeCollected.get(userWalk).size(); i++) {
//					System.out.println(i);
					// viene valutata l'affermazione del probe e i valori di FA e FR vengono aggiornati di conseguenza
					checkIdentityClaim(userWalk, i, threshold);
				}
			}
			
			// per ogni threshold, vengono calcolati FAR e FRR
			// i valori utilizzati per il calcolo dei rate, sono dati rispettivamente dal numero totale di dichiarazioni false e genuine di tutti i probe
			far.add(localFA/14112.0);
			frr.add(localFR/294.0);
			
			totalGA += localGA;
			totalFA += localFA;
			totalGR += localGR;
			totalFR += localFR;
		}
		
		// vengono calcolati i valori per ROC
		roc = frr.stream().map(e -> 1 - e).collect(Collectors.toList());
	}
	
	/**
	 * esegue la verifica sull'identità del probe in base al threshold
	 * @param userWalk è la camminata del probe e contiene informazioni circa la sua reale identità 
	 * @param identity è l'identità dichiarata dal probe
	 * @param threshold è il valore fissato per l'accettazione del subject
	 */
	private void checkIdentityClaim(UserWalk userWalk, int identity, double threshold) {
		
		// si controlla che il valore della distanza tra il probe e il subject della gallery sia minore del threshold
		// identity è un valore originato da un ciclo con indice iniziale 1, ma la lista dei valori delle distanze per i subject hanno indice iniziale 0
		// è quindi necessario decrementare localmente il valore dell'identità per avere un confronto reale
		if(probeCollected.get(userWalk).get(identity-1) <= threshold) {
			// viene valutata la genuinità della dichiarazione sull'identità
			// identity è un valore originato da un ciclo con indice iniziale 0, ma i subject sono identificati da un intero con indice iniziale 1
			// è quindi necessario incrementare il valore dell'identità per avere un confronto reale
			if(userWalk.getUser() != identity) {
				// il probe è un impostore e si ha un falso positivo
				localFA++;
			} else {
				localGA++;
			}
		} else {
			// la distanza misurata è inferiore al threshold
			if(userWalk.getUser() == identity) {
				// il probe è genuino e si ha un falso negativo
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
