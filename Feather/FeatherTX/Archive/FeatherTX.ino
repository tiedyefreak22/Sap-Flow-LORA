#include <SPI.h>
#include <RH_RF95.h>
#include <RHReliableDatagram.h>

#define CLIENT_ADDRESS 1
#define SERVER_ADDRESS 2

#define RFM95_CS    8
#define RFM95_INT   3
#define RFM95_RST   4

// Change to 434.0 or other frequency, must match RX's freq!
#define RF95_FREQ 915.0
int led = LED_BUILTIN;
int analogPin = A1; // sensor signal pin
int val = 0;  // variable to store the value read

// Singleton instance of the radio driver
RH_RF95 rf95(RFM95_CS, RFM95_INT);

// Class to manage message delivery and receipt, using the driver declared above
RHReliableDatagram manager(rf95, CLIENT_ADDRESS);

void setup() {
  pinMode(RFM95_RST, OUTPUT);
  digitalWrite(RFM95_RST, HIGH);
  
  // set LED to be an output pin
  pinMode(led, OUTPUT);

  Serial.begin(115200);
  while (!Serial) delay(1);
  delay(100);

  Serial.println("Feather LoRa TX Test!");

  // manual reset
  digitalWrite(RFM95_RST, LOW);
  delay(10);
  digitalWrite(RFM95_RST, HIGH);
  delay(10);

  while (!rf95.init()) {
    Serial.println("LoRa radio init failed");
    Serial.println("Uncomment '#define SERIAL_DEBUG' in RH_RF95.cpp for detailed debug info");
    while (1);
  }
  Serial.println("LoRa radio init OK!");

  // Defaults after init are 434.0MHz, modulation GFSK_Rb250Fd250, +13dbM
  if (!rf95.setFrequency(RF95_FREQ)) {
    Serial.println("setFrequency failed");
    while (1);
  }
  Serial.print("Set Freq to: "); Serial.println(RF95_FREQ);

  // Defaults after init are 434.0MHz, 13dBm, Bw = 125 kHz, Cr = 4/5, Sf = 128chips/symbol, CRC on

  // The default transmitter power is 13dBm, using PA_BOOST.
  // If you are using RFM95/96/97/98 modules which uses the PA_BOOST transmitter pin, then
  // you can set transmitter powers from 5 to 23 dBm:
  rf95.setTxPower(23, false);
}

int16_t packetnum = 0;  // packet counter, we increment per xmission

void loop() {
  delay(10000); // Wait 1 second between transmits, could also 'sleep' here!
  Serial.println("Transmitting..."); // Send a message to rf95_server
  digitalWrite(led, HIGH);    // turn the LED on

  val = analogRead(analogPin);
  char radiopacket[20];
  sprintf(radiopacket, "%d", val);
  itoa(packetnum++, radiopacket+13, 10);
  Serial.print("Sending "); Serial.println(radiopacket);
  radiopacket[19] = 0;

  Serial.println("Sending...");
  delay(10);


  if (manager.sendtoWait(buf, strlen((uint8_t *)radiopacket, SERVER_ADDRESS))
	{
		// Now wait for a reply from the server
		uint8_t len = sizeof(buf);
		uint8_t from;
		if (manager.recvfromAckTimeout(buf, &len, 2000, &from))
		{
			buf[len] = 0;
			Serial.printf("got reply from 0x%02x rssi=%d %s", from, driver.lastRssi(), (char *) buf);
		}
		else
		{
			Serial.println("No reply, is rf95_reliable_datagram_server running?");
		}
	}
	else
		Serial.println("sendtoWait failed");
  digitalWrite(led, LOW);    // turn the LED off
	delay(500);
}
