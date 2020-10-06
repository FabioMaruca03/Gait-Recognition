package ImageAnalysis;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import javax.imageio.ImageIO;

import DataCenter.PersonWalk;

public class ImageCreator {

	/**
	 * crea l'immagine e la salva
	 * 
	 * @param width e' la lunghezza che avra' l'immagine
	 * @param source e' il path contenente il dataset parallelo
	 * @param destination e' il path in cui andare a salvare l'immagine
	 * @param userWalk contiene i dati della camminata del subject
	 */
	public void createImage(Double width, String source, String destination, UserWalk userWalk) {

		// legge il file del dataset
		try (BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(source + "\\imageData.csv")))) {
			// istanzia l'oggetto BufferedImage che conterra' la nuova immagine
			// i parametri di input sono rispettivamente la lunghezza, l'altezza e il ColorSpace utilizzato
			// la lunghezza dei passi e' fissa, mentre il numero degli stessi e' preso dinamicamente per ogni immagine da userWalk
			// stepsList e' una lista di oggetti custom contenuta in userWalk con le informazioni sui passi della camminata
			// il numero di elementi della lista rappresenta il numero di passi che compongono la camminata stessa 
			BufferedImage image = new BufferedImage(width.intValue(), userWalk.getStepsList().size(), BufferedImage.TYPE_INT_RGB);

			String currentLine;
			int x = 0;
			int y = 0;

			// scorre il file di input
			while ((currentLine = bufferedReader.readLine()) != null) {
				String[] line = currentLine.split(";");

				// ogni riga dell'immagine viene riempita in base ai dati del singolo passo
				if (x < width) {
					// le informazioni dei canali vengono usate per creare l'oggetto Color, che rappresenta il pixel
					Color color = new Color(Integer.parseInt(line[1]), Integer.parseInt(line[2]), Integer.parseInt(line[3]));
					// viene salvato il pixel nella posizione indicata
					image.setRGB(x, y, color.getRGB());
					x++;
				} else {
					// una volta completata la riga si passa alla successiva
					y++;
					Color color = new Color(Integer.parseInt(line[1]), Integer.parseInt(line[2]), Integer.parseInt(line[3]));
					image.setRGB(0, y, color.getRGB());
					x = 1;
				}
			}

			// viene creata la cartella contenente le immagini
			File directory = new File(destination);
			if (!directory.exists()) {
				directory.mkdirs();
			}

			// l'immagine viene salvata in un nuovo file
			File imageFile = new File(destination + userWalk.getWalk() + ".png");
			ImageIO.write(image, "png", imageFile);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	public void createImage2(double width, String destination, PersonWalk personWalk, int walkId) {


		// istanzia l'oggetto BufferedImage che conterra' la nuova immagine
		// i parametri di input sono rispettivamente la lunghezza, l'altezza e il ColorSpace utilizzato
		// la lunghezza dei passi e' fissa, mentre il numero degli stessi e' preso dinamicamente per ogni immagine da userWalk
		// stepsList e' una lista di oggetti custom contenuta in userWalk con le informazioni sui passi della camminata
		// il numero di elementi della lista rappresenta il numero di passi che compongono la camminata stessa 
		BufferedImage image = new BufferedImage((int)width, personWalk.getSteps().size(), BufferedImage.TYPE_INT_RGB);

		String currentLine;
		int x = 0;
		int y = 0;

		for (int i=0;i<personWalk.getArrayDoubleInterpolationX().size();i++) {
			// ogni riga dell'immagine viene riempita in base ai dati del singolo passo
			if (x < width) {
				// le informazioni dei canali vengono usate per creare l'oggetto Color, che rappresenta il pixel
				Color color = new Color(personWalk.getArrayDoubleInterpolationX().get(i).intValue(),personWalk.getArrayDoubleInterpolationY().get(i).intValue(),personWalk.getArrayDoubleInterpolationZ().get(i).intValue());
				// viene salvato il pixel nella posizione indicata
				image.setRGB(x, y, color.getRGB());
				x++;
			} else {
				// una volta completata la riga si passa alla successiva
				y++;
				Color color = new Color(personWalk.getArrayDoubleInterpolationX().get(i).intValue(),personWalk.getArrayDoubleInterpolationY().get(i).intValue(),personWalk.getArrayDoubleInterpolationZ().get(i).intValue());
				image.setRGB(0, y, color.getRGB());
				x = 1;
			}
		}
		// viene creata la cartella contenente le immagini
		File directory = new File(destination);
		if (!directory.exists()) {
			directory.mkdirs();
		}

		// l'immagine viene salvata in un nuovo file
		File imageFile = new File(destination + walkId + ".png");
		try {
			ImageIO.write(image, "png", imageFile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}



}
