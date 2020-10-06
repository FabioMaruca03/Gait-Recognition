package ImageAnalysis;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class LBP {
	
	/**
	 * genera nuove immagini a partire da quelle originali andando ad applicare LBP
	 * @param source e' il path del file contenente l'immagine a cui applicare LBP. Lo stesso source viene usato per salvare la nuova immagine con il nome modificato in modo da evitare sovrascritture
	 * @param userWalk contiente informazioni sulla camminata
	 * @param width e' la lunghezza del passo
	 */
	public void calculateLBP(String source, UserWalk userWalk, Double width) {
		
		try {
			// istanzia BufferedImage per andare a creare la nuova immagine per LBP con la stessa struttura dell'immagine di partenza
			// viene caricata l'immagine originale
			BufferedImage lbpImage = new BufferedImage(width.intValue(), userWalk.getStepsList().size(), BufferedImage.TYPE_INT_RGB);
			BufferedImage image = ImageIO.read(new File(source + ".png"));
			
			// applica lbp ad ogni pixel
			for(int y = 0; y < image.getHeight(); y++) {
				for(int x = 0; x < image.getWidth(); x++) {
					
					// queste stringhe conterranno i risultati di LBP applicato al pixel in questione
					String lbpRed = "";
					String lbpGreen = "";
					String lbpBlue = "";
					
					// viene caricato il pixel dall'immagine
					int pixel = image.getRGB(x, y);
					
					// le informazioni del pixel vengono scomposte nei canali R, G e B
					int r = (pixel>>16) & 0xff;
					int g = (pixel>>8) & 0xff;
					int b = pixel & 0xff;
					
					// vengono registrati i valori ottenuti dall'analisi di ogni pixel nell'intorno del pixel di partenza
					// per mantenere una struttura circolare nel salvataggio dei valori, la funzione viene chiamata più volte passando parametri diversi per l'offset di x e y
					// l'offset va da -1 a 1 e si esclude la coppia [0, 0], in modo da non confrontare il pixel con se' stesso
					String[] resultUpLeft = generateStringLBPH(image, x, y, -1, -1, lbpRed, lbpGreen, lbpBlue, r, g, b);
					lbpRed = resultUpLeft[0];
					lbpGreen = resultUpLeft[1];
					lbpBlue = resultUpLeft[2];
					
					String[] resultUpCenter = generateStringLBPH(image, x, y, 0, -1, lbpRed, lbpGreen, lbpBlue, r, g, b);
					lbpRed = resultUpCenter[0];
					lbpGreen = resultUpCenter[1];
					lbpBlue = resultUpCenter[2];
					
					String[] resultUpRight = generateStringLBPH(image, x, y, 1, -1, lbpRed, lbpGreen, lbpBlue, r, g, b);
					lbpRed = resultUpRight[0];
					lbpGreen = resultUpRight[1];
					lbpBlue = resultUpRight[2];
					
					String[] resultCenterLeft = generateStringLBPH(image, x, y, -1, 0, lbpRed, lbpGreen, lbpBlue, r, g, b);
					lbpRed = resultCenterLeft[0];
					lbpGreen = resultCenterLeft[1];
					lbpBlue = resultCenterLeft[2];
					
					String[] resultCenterRight = generateStringLBPH(image, x, y, 1, 0, lbpRed, lbpGreen, lbpBlue, r, g, b);
					lbpRed = resultCenterRight[0];
					lbpGreen = resultCenterRight[1];
					lbpBlue = resultCenterRight[2];
					
					String[] resultDownLeft = generateStringLBPH(image, x, y, -1, 1, lbpRed, lbpGreen, lbpBlue, r, g, b);
					lbpRed = resultDownLeft[0];
					lbpGreen = resultDownLeft[1];
					lbpBlue = resultDownLeft[2];
					
					String[] resultDownCenter = generateStringLBPH(image, x, y, 0, 1, lbpRed, lbpGreen, lbpBlue, r, g, b);
					lbpRed = resultDownCenter[0];
					lbpGreen = resultDownCenter[1];
					lbpBlue = resultDownCenter[2];
					
					String[] resultDownRight = generateStringLBPH(image, x, y, 1, 1, lbpRed, lbpGreen, lbpBlue, r, g, b);
					lbpRed = resultDownRight[0];
					lbpGreen = resultDownRight[1];
					lbpBlue = resultDownRight[2];
					
					// le informazioni per il colore vengono prese a partire dall'array di valori e riportandoli in interi in base 2
					Color color = new Color(Integer.parseInt(lbpRed, 2), Integer.parseInt(lbpGreen, 2), Integer.parseInt(lbpBlue, 2));
					lbpImage.setRGB(x, y, color.getRGB());
				}
			}
			
			// viene creata l'immagine
			File imageFile = new File(source + "lbp.png");
			ImageIO.write(lbpImage, "png", imageFile);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * confronta il pixel con il vicino selezionato e ne restituisce la relazione
	 * @param image e' l'immagine caricata
	 * @param x e' la coordinata x del pixel in esame
	 * @param y e' la coordinata y del pixel in esame
	 * @param offsetX e' il modificatore applicata alla coordinata x del pixel per ottenere la coordinata del pixel vicino
	 * @param offsetY e' il modificatore applicata alla coordinata x del pixel per ottenere la coordinata del pixel vicino
	 * @param lbpRed e' la stringa contenente il risultato dei confronti gia' eseguiti per il canale R
	 * @param lbpGreen e' la stringa contenente il risultato dei confronti gia' eseguiti per il canale G
	 * @param lbpBlue e' la stringa contenente il risultato dei confronti gia' eseguiti per il canale B
	 * @param r e' il valore del canale R del pixel in esame
	 * @param g e' il valore del canale G del pixel in esame
	 * @param b e' il valore del canale B del pixel in esame
	 * @return l'insieme delle stringhe lbpRed, lbpGreen e lbpBlue con l'aggiunta del risultato dell'iterazione corrente
	 */
	private String[] generateStringLBPH(BufferedImage image, int x, int y, int offsetX, int offsetY, String lbpRed, String lbpGreen, String lbpBlue, int r, int g, int b) {
		
		// si controlla che il pixel in esame non sia al bordo dell'immagine
		if(x + offsetX >= 0 && y + offsetY >= 0 && x + offsetX < image.getWidth() && y + offsetY < image.getHeight()) {
			
			// viene estratto il pixel del vicino
			int lbpPixel = image.getRGB(x + offsetX, y + offsetY);
			
			// il pixel viene scomposto
			int lbpR = (lbpPixel>>16) & 0xff;
			int lbpG = (lbpPixel>>8) & 0xff;
			int lbpB = lbpPixel & 0xff;
			
			// si confronta il valore per ogni canale del vicino con quello corrispondente del pixel originale
			// la strigna risultante viene modificata di conseguenza
			if(lbpR >= r) {
				lbpRed += "1";
			} else {
				lbpRed += "0";
			}
			
			if(lbpG >= g) {
				lbpGreen += "1";
			} else {
				lbpGreen += "0";
			}
			
			if(lbpB >= b) {
				lbpBlue += "1";
			} else {
				lbpBlue += "0";
			}
		} else {
			// se il pixel originale appartiene al bordo dell'immagine, viene assegnato di default il valore 0
			lbpRed += "0";
			lbpGreen += "0";
			lbpBlue += "0";
		}
		
		String[] result = {lbpRed, lbpGreen, lbpBlue};
		
		return result;
	}
}
