package ImageAnalysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

public class DataPreparation {

	int allStepsLength = 0;
	List<Double> minMax = new ArrayList<>();
	Double stdevSteps = 0.0;

	Double meanStepLength = 0.0;
	public Double getMeanStepLength() { return meanStepLength; }
	
	double maxValueX = 0;
	double minValueX = 0;
	double maxValueY = 0;
	double minValueY = 0;
	double maxValueZ = 0;
	double minValueZ = 0;
	
	/**
	 * preapara i dati per creare le immagini
	 * @return la struttura delle camminate
	 */
	public List<UserWalk> getUsersWalks() {

		List<UserWalk> usersWalks = new ArrayList<>();
		Double sumMeanSteps = 0.0;

		// recupera i file
		for (int user = 1; user <= 153; user++) {
			for (int walk = 0; walk < 12; walk++) {

				String subj = "subj_";

				if (user / 100 >= 1) {
					subj += String.valueOf(user);
				} else if (user / 10 >= 1) {
					subj += "0" + user;
				} else {
					subj += "00" + user;
				}

				FileReader fileReader;
				FileReader fileReaderStepsPreparation;

				try {
					// crea l'oggetto userWalk e setta subject e camminata
					UserWalk userWalk = new UserWalk();
					userWalk.setUser(user);
					userWalk.setWalk(walk);

					String currentLine;

					// inizializza le liste che andranno salvate in userWalk
					List<Integer> walkStartIndex = new ArrayList<>();
					List<Integer> walkEndIndex = new ArrayList<>();
					List<Step> stepList = new ArrayList<>();

					// necessario per accedere al file csv corretto
					String trimmeredWalk = String.valueOf(walk);
					if (trimmeredWalk.length() > 1) {
						trimmeredWalk = trimmeredWalk.substring(1);
					}

					// vengono caricati i file necessari a processare i dati
					// fileReader: legge gli indici dei passi per determinare inizio e fine di ognuno
					// fileReaderStepsPreparation: legge il file con i valori per associare ogni step alla camminata

					// pc casa
					fileReader = new FileReader(new File(Constants.datasetSource +File.pathSeparator+ subj +File.pathSeparator+ "walk" + walk  +File.pathSeparator+"data"+File.pathSeparator+"stepsIndex.txt"));
					fileReaderStepsPreparation = new FileReader(new File(Constants.datasetSource + subj +File.pathSeparator+ "walk" + walk +File.pathSeparator+ subj + trimmeredWalk + ".csv"));

					// pc lavoro
					// fileReader = new FileReader(new File("C:\\Users\\ggiorgi\\Desktop\\Tirocinio\\DatasetChinaComplete-Segmented\\" + subj + "\\walk" + walk + "\\data\\stepsIndex.txt"));
					// fileReaderStepsPreparation = new FileReader(new File("C:\\Users\\ggiorgi\\Desktop\\Tirocinio\\DatasetChinaComplete-Segmented\\" + subj + "\\walk" + walk + "\\" + subj + trimmeredWalk + ".csv"));

					// associa gli indici dei passi alla camminata
					try (BufferedReader bufferedReader = new BufferedReader(fileReader)) {
						while ((currentLine = bufferedReader.readLine()) != null) {
							String[] line = currentLine.split(";");
							// aggiunge l'indice iniziale e finale di ogni passo
							walkStartIndex.add(Integer.valueOf(line[0]));
							walkEndIndex.add(Integer.valueOf(line[1]));
						}

						// salva la lista degli indici in userWalk
						userWalk.setStepsStartIndex(walkStartIndex);
						userWalk.setEndIndex(walkEndIndex);

						// estrae i singoli step e li associa alle camminate, insieme alla lunghezza di ognuno
						try (BufferedReader bufferedReaderStepsPreparation = new BufferedReader(fileReaderStepsPreparation)) {

							List<String> readFile = new ArrayList<>();
							while ((currentLine = bufferedReaderStepsPreparation.readLine()) != null) {
								// viene copiato il contenuto del file in una lista in modo da poter leggere pi� volte il file senza dover aprire ogni volta il buffer
								readFile.add(currentLine);
							}

							// crea uno Step per ogni passo della camminata
							for (int i = 0; i < userWalk.getStepsStartIndex().size(); i++) {
								Step step = new Step();
								List<String> stepString = new ArrayList<>();

								for (String readLine : readFile) {
									String[] line = readLine.split(";");

									// se il passo � all'interno di stepIndex viene aggiunto, altrimenti viene scartato
									if (userWalk.getStepsStartIndex().get(i) <= Integer.valueOf(line[0]) && Integer.valueOf(line[0]) <= userWalk.getEndIndex().get(i)) {
										String result = line[0] + ";" + Double.valueOf(line[1]) + ";" + Double.valueOf(line[2]) + ";" + Double.valueOf(line[3]) + ";";
										stepString.add(result);
									}
								}

								// l'ultimo valore della lunghezza del passo � esclusivo, quindi viene scartato
								stepString.remove(stepString.size() - 1);
								// viene salvato il passo e la lunghezza
								step.setStep(stepString);
								step.setLength(stepString.size());
								// il passo viene aggiunto alla lista dei passi della camminata
								stepList.add(step);
							}

						} catch (IOException e) {
							e.printStackTrace();
						}

						// vengono aggiunti i passi alla camminata
						userWalk.setStepsList(stepList);
						// la camminata viene aggiunta alla lista completa
						usersWalks.add(userWalk);

					} catch (IOException e) {
						e.printStackTrace();
					}
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}
			}
		}

		// calcola la lunghezza totale dei passi
		for (UserWalk walk : usersWalks) {
			for (Step step : walk.getStepsList()) {
				allStepsLength += step.getLength();
			}
		}

		Double totalSteps = 0.0;
		
		// calcola il numero totale di passi tra tutte le camminate
		for(UserWalk userWalk: usersWalks) {
			totalSteps += userWalk.getStepsList().size();
		}
		
		// calcola la lunghezza media di ogni passo
		meanStepLength = allStepsLength / totalSteps;
		
		// calcola la deviazione standard
		for(UserWalk userWalk : usersWalks) {
			for(Step step : userWalk.getStepsList()) {
				sumMeanSteps += Math.pow(step.getLength()/meanStepLength, 2);
			}
		}
		
		stdevSteps = Math.sqrt(sumMeanSteps / totalSteps);
		
		return usersWalks;
	}

	/**
	 * calcola i valori interpolati a partire dal dataset originale
	 * @param folderPath e' il path da cui andare a creare il file con i dati interpolati  
	 * @param userWalk contiene i dati dei passi da interpolare
	 */
	public void dataInterpolation(String folderPath, UserWalk userWalk) {

		// istanzia l'outputStream per scrivere il file 
		try (FileOutputStream fileWriter = new FileOutputStream(folderPath + "\\imageDataInterpolated.csv")) {
			// itera su tutti gli step che compongono la camminata
			for (Step step : userWalk.getStepsList()) {
				
				// il rapporto tra la lunghezza del passo preso in considerazione e il valore medio verr� usato per andare a determinare i nuovi valori interpolati
				Double ratio = Double.valueOf(step.getLength() - 1) / meanStepLength;

				// vengono istanziati degli array che conterranno i valori che compongono il passo 
				double[] indexArray = new double[step.getLength()];
				double[] xArray = new double[step.getLength()];
				double[] yArray = new double[step.getLength()];
				double[] zArray = new double[step.getLength()];

				// vengono riempiti gli array, ognuno con il valore di un asse
				for(int i = 0; i < step.getLength(); i++) {
					String[] line = step.getStep().get(i).split(";");
					xArray[i] = Double.valueOf(line[1]);
					yArray[i] = Double.valueOf(line[2]);
					zArray[i] = Double.valueOf(line[3]);
				}

				// indexArray viene riempito assegnando ad ogni indice il proprio valore numerico
				for (int i = 0; i < indexArray.length; i++) {
					indexArray[i] = Double.valueOf(i);
				}

				// istanzia l'interpolatore lineare
				LinearInterpolator li = new LinearInterpolator();

				// calcola le funzioni di interpolazione
				PolynomialSplineFunction xPsf = li.interpolate(indexArray, xArray);
				PolynomialSplineFunction yPsf = li.interpolate(indexArray, yArray);
				PolynomialSplineFunction zPsf = li.interpolate(indexArray, zArray);

				// calcola i valori nei singoli punti interpolati e li scrive sul file
				for (int i = 0; i < meanStepLength; i++) {
					String result = (i * ratio) + ";" + xPsf.value(i * ratio) + ";" + yPsf.value(i * ratio) + ";" + zPsf.value(i * ratio) + "\n";
					fileWriter.write(result.getBytes());
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
		
	/**
	 * crea un dataset parallelo a quello iniziale
	 * @param source e' il path da cui caricare il file con i dati interpolati
	 */
	public void createImageDataset(String source) {
		
		// apre il file con i dati interpolati
		try(BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(source + "\\imageDataInterpolated.csv")));
				// istanzia l'OutputStream per salvare i dati normalizzati
				FileOutputStream fileWriter = new FileOutputStream(source + "\\imageData.csv");) {

			String currentLine;

			// itera sul file aperto
			while ((currentLine = bufferedReader.readLine()) != null) {
				String[] line = currentLine.split(";");

				// calcola i valori normalizzati, li porta in una scala [0, 255] e li salva
				// minMax e' una lista i cui indici pari contengono i massimi per ogni asse, gli indici dispari contengono i minimi
				int x = (int) Math.round(((Double.parseDouble(line[1]) - minMax.get(1)) / (minMax.get(0) - minMax.get(1))) * 255);
				int y = (int) Math.round(((Double.parseDouble(line[2]) - minMax.get(3)) / (minMax.get(2) - minMax.get(3))) * 255);
				int z = (int) Math.round(((Double.parseDouble(line[3]) - minMax.get(5)) / (minMax.get(4) - minMax.get(5))) * 255);

				String result = line[0] + ";" + x + ";" + y + ";" + z + ";" + (x + y + z) / 3 + "\n";
				fileWriter.write(result.getBytes());
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	// calcola minimo e massimo
	public void generateMinMax(List<UserWalk> usersWalks) {

		for (UserWalk userWalk : usersWalks) {

			// crea la stringa per accedere alla folder dell'utente
			String subj = "subj_";

			if (userWalk.getUser() / 100 >= 1) {
				subj += String.valueOf(userWalk.getUser());
			} else if (userWalk.getUser() / 10 >= 1) {
				subj += "0" + userWalk.getUser();
			} else {
				subj += "00" + userWalk.getUser();
			}

			// apre il file con i dati interpolati
			try (BufferedReader bufferedReader = new BufferedReader(
					new FileReader(new File("C:\\Users\\kenam\\Desktop\\Tirocinio\\DatasetChinaComplete-Segmented\\"
							+ subj + "\\walk" + userWalk.getWalk() + "\\imageDataInterpolated.csv")))) {

				String currentLine;

				// legge le righe
				while ((currentLine = bufferedReader.readLine()) != null) {
					String[] line = currentLine.split(";");

					// aggiorna minimo e massimo per ogni asse
					maxValueX = Math.max(maxValueX, Double.parseDouble(line[1]));
					minValueX = Math.min(minValueX, Double.parseDouble(line[1]));
					maxValueY = Math.max(maxValueY, Double.parseDouble(line[2]));
					minValueY = Math.min(minValueY, Double.parseDouble(line[2]));
					maxValueZ = Math.max(maxValueZ, Double.parseDouble(line[3]));
					minValueZ = Math.min(minValueZ, Double.parseDouble(line[3]));
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		// setta minimo e massimo in variabili globali
		minMax.add(maxValueX);
		minMax.add(minValueX);
		minMax.add(maxValueY);
		minMax.add(minValueY);
		minMax.add(maxValueZ);
		minMax.add(minValueZ);
	}
}