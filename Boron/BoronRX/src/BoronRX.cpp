/* 
 * Project myProject
 * Author: Your Name
 * Date: 
 * For comprehensive documentation and examples, please visit:
 * https://docs.particle.io/firmware/best-practices/firmware-template/
 */

// Include Particle Device OS APIs
#include "Particle.h"
#include <SPI.h>
#include <RH_RF95.h>

// Let Device OS manage the connection to the Particle Cloud
SYSTEM_MODE(AUTOMATIC);

// Run the application and system concurrently in separate threads
SYSTEM_THREAD(ENABLED);

// Show system, cloud connectivity, and application logs over USB
// View logs with CLI using 'particle serial monitor --follow'
// SerialLogHandler logHandler(LOG_LEVEL_INFO);

#define RFM95_CS A5
#define RFM95_RST D5
#define RFM95_INT D2
#define RF95_FREQ 915.0
char* temp_feather_1;

// Singleton instance of the radio driver
RH_RF95 rf95(RFM95_CS, RFM95_INT);

// setup() runs once, when the device is first turned on
void setup() {
  pinMode(RFM95_RST, OUTPUT);
  digitalWrite(RFM95_RST, HIGH);
  //Particle.variable("temp_feather_1", &temp_feather_1, DOUBLE);
  Particle.connect();

  Serial.begin(115200);
  delay(2000);

  Serial.println("Feather LoRa RX Test!");

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

  // Defaults after init are 434.0MHz, modulation GFSK_Rb250Fd250, +13dbM
  if (!rf95.setFrequency(RF95_FREQ)) {
    Serial.println("setFrequency failed");
    while (1);
  }
  Serial.print("Set Freq to: "); Serial.println(RF95_FREQ);
  rf95.setTxPower(23, false);
}

void loop() {
  if (rf95.available()) {
    // Should be a message for us now
    uint8_t buf[RH_RF95_MAX_MESSAGE_LEN];
    uint8_t len = sizeof(buf);

    if (rf95.recv(buf, &len)) {
      temp_feather_1 = (char*)buf;
      //RH_RF95::printBuffer("Received: ", temp_feather_1, len);
      Serial.print("Got: ");
      Serial.println(temp_feather_1);
      // if (Particle.connected()) {
      //   bool success = Particle.publish("temp_f_1", temp_feather_1, WITH_ACK);
      //   if (success) {
      //     Serial.println("Published data to cloud");
      //   }
      // }
      // else {
      //   Serial.println("Could not publish data to cloud");
      // }
      Serial.print("RSSI: ");
      Serial.println(rf95.lastRssi(), DEC);

      // Send a reply
        char radiopacket[20] = {' '};
        sprintf(radiopacket, temp_feather_1);
        rf95.send((uint8_t *)radiopacket, strlen(radiopacket)+1);
        delay(10);
        rf95.waitPacketSent();
        Serial.println("Sent a reply");
    } else {
        Serial.println("Receive failed");
    }
  }
}
