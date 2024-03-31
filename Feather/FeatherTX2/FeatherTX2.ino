// Feather9x_TX
// -*- mode: C++ -*-
// Example sketch showing how to create a simple messaging client (transmitter)
// with the RH_RF95 class. RH_RF95 class does not provide for addressing or
// reliability, so you should only use RH_RF95 if you do not need the higher
// level messaging abilities.
// It is designed to work with the other example Feather9x_RX

#include <SPI.h>
#include <RH_RF95.h>
#include <RHReliableDatagram.h>

#define RFM95_CS    8
#define RFM95_INT   3
#define RFM95_RST   4

#define CLIENT_ADDRESS 1
#define SERVER_ADDRESS 2

int led = LED_BUILTIN;
int sense1 = A0; // sensor signal pin1
int sense2 = A2; // sensor signal pin
float val1 = 0.0;  // variable to store the value read
float val2 = 0.0;  // variable to store the value read

// Change to 434.0 or other frequency, must match RX's freq!
#define RF95_FREQ 915.0

// Singleton instance of the radio driver
RH_RF95 rf95(RFM95_CS, RFM95_INT);
RHReliableDatagram manager(rf95, CLIENT_ADDRESS);

void setup() {
  pinMode(led, OUTPUT);
  pinMode(RFM95_RST, OUTPUT);
  digitalWrite(RFM95_RST, HIGH);

  Serial.begin(115200);
  delay(2000);

  Serial.println("Feather LoRa TX Test!");

  // manual reset
  digitalWrite(RFM95_RST, LOW);
  delay(100);
  digitalWrite(RFM95_RST, HIGH);
  delay(10);

  while (!rf95.init()) {
    Serial.println("LoRa radio init failed");
    Serial.println("Uncomment '#define SERIAL_DEBUG' in RH_RF95.cpp for detailed debug info");
    while (1);
  }
  Serial.println("LoRa radio init OK!");

  if (!manager.init())
		Serial.println("init failed");
  Serial.println("Manager init OK!");

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

int16_t counter = 0;  // packet counter, we increment per xmission
uint8_t buf[RH_RF95_MAX_MESSAGE_LEN];

void loop() {
  delay(5000);
  
  val1 = analogRead(sense1);

  // Gymnastics required to get analogRead val into float
  char radiopacket[20] = {' '};
  int val1_int = (int) val1;
  float val1_float = (abs(val1) - abs(val1_int)) * 100000;
  int val1_fra = (int)val1_float;
  sprintf(radiopacket, "%d.%d", val1_int, val1_fra);
  
  itoa(counter++, radiopacket+13, 10);
  Serial.print("Sending "); Serial.println(radiopacket);
  digitalWrite(led, HIGH);    // turn the LED on

  delay(10);
  Serial.println("Sending to rf95_reliable_datagram_server");

	// Send a message to manager_server
	snprintf((char *)buf, sizeof(buf), "to server counter=%d", ++counter);

	if (manager.sendtoWait(buf, strlen((char *)buf), SERVER_ADDRESS))
	{
		// Now wait for a reply from the server
		uint8_t len = sizeof(buf);
		uint8_t from;
		if (manager.recvfromAckTimeout(buf, &len, 2000, &from))
		{
			buf[len] = 0;
			Serial.printf("got reply from 0x%02x rssi=%d %s", from, rf95.lastRssi(), (char *) buf);
		}
		else
		{
			Serial.println("No reply, is rf95_reliable_datagram_server running?");
		}
	}
	else
		Serial.println("sendtoWait failed");
	delay(500);

  digitalWrite(led, LOW);    // turn the LED off
}