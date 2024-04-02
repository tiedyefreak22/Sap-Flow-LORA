#include "Particle.h"
#include <SPI.h>
#include <RH_RF95.h>
#include <RHReliableDatagram.h>
#define CLIENT_ADDRESS 1
#define SERVER_ADDRESS 2

// Let Device OS manage the connection to the Particle Cloud
SYSTEM_MODE(AUTOMATIC);

// Run the application and system concurrently in separate threads
SYSTEM_THREAD(ENABLED);

SerialLogHandler logHandler;

// How often to publish a value
const std::chrono::milliseconds publishPeriod = 30s;

// The event name to publish with
const char *eventName = "sheetTest1";

unsigned long lastPublish;
int counter = 0;

void publishTest();

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
RHReliableDatagram manager(rf95, SERVER_ADDRESS);

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

  if (!manager.init())
		Serial.println("Manager init failed");
  Serial.println("Manager init OK!");

  // Defaults after init are 434.0MHz, modulation GFSK_Rb250Fd250, +13dbM
  if (!rf95.setFrequency(RF95_FREQ)) {
    Serial.println("setFrequency failed");
    while (1);
  }
  Serial.print("Set Freq to: "); Serial.println(RF95_FREQ);
  rf95.setTxPower(23, false);
}

uint8_t buf[RH_RF95_MAX_MESSAGE_LEN];

void loop() {
//   if (rf95.available()) {
  if (manager.available()) {
    //if (rf95.recv(buf, &len)) {
    temp_feather_1 = (char*)buf;
    // Wait for a message addressed to us from the client
    uint8_t len = sizeof(buf);
    uint8_t from;
    if (manager.recvfromAck(buf, &len, &from))
    {
        buf[len] = 0;
        Serial.printlnf("got packet from 0x%02x rssi=%d %s", from, rf95.lastRssi(), temp_feather_1);

        int request = 0;
        char *cp = strchr(temp_feather_1, '=');
        if (cp) {
            request = atoi(cp + 1);
        }

        snprintf(temp_feather_1, sizeof(buf), "request=%d rssi=%d", request, rf95.lastRssi());

        // Send a reply back to the originator client
        if (!manager.sendtoWait(buf, strlen(temp_feather_1), from))
            Serial.println("sendtoWait failed");

        if (Particle.connected())
        {
          if (millis() - lastPublish >= publishPeriod.count()) {
              lastPublish = millis();

              publishTest();
          }
        }
    } else {
        Serial.println("Receive failed");
    }
  }
}

void publishTest() {
    char buf[128];

    snprintf(buf, sizeof(buf), "[%d,%d]", ++counter, rand());

    Particle.publish(eventName, buf, PRIVATE);
    Log.info("published: %s", buf);
}