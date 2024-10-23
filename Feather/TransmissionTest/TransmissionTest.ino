// f1 MAC: 98:76:B6:12:99:a5
#include <RHReliableDatagram.h>
#include <RH_RF95.h>
#include <SPI.h>
#include <SD.h>

#define RFM95_CS    8
#define RFM95_INT   3
#define RFM95_RST   4
#define LED_RED     LED_BUILTIN
#define SD_CS       10

#define CLIENT_ADDRESS 1
#define SERVER_ADDRESS 2

#define SYNC_WORD 0x34

// Should these be switched to #define?
int led = LED_BUILTIN;
int sense1 = A0; // sensor signal pin 1
int sense2 = A2; // sensor signal pin 2
float val1 = 0.0;  // variable to store the value read
float val2 = 0.0;  // variable to store the value read
int REPORT_INT = 10.0; // interval at which to record and report data

// Change to 434.0 or other frequency, must match RX's freq!
#define RF95_FREQ 915.0

// Singleton instance of the radio driver
RH_RF95 rf95(RFM95_CS, RFM95_INT);
RHReliableDatagram manager(rf95, CLIENT_ADDRESS);

File logfile;

void setup() {
  pinMode(led, OUTPUT);
  pinMode(RFM95_RST, OUTPUT);

  // RF95 Initialization
  // manual reset
  digitalWrite(RFM95_RST, HIGH);
  digitalWrite(RFM95_RST, LOW);
  delay(100);
  digitalWrite(RFM95_RST, HIGH);
  delay(10);

  Serial.begin(9600);
  delay(2000);

  // if (!manager.init())
	// 	Serial.println("Manager init failed");
  //   error(2);
  //   while(1);
  // Serial.println("Manager init OK!");

  while (!rf95.init()) {
    Serial.println("LoRa radio init failed");
    Serial.println("Uncomment '#define SERIAL_DEBUG' in RH_RF95.cpp for detailed debug info");
    error(1);
    while (1);
  }
  Serial.println("LoRa radio init OK!");

  // Defaults after init are 434.0MHz, modulation GFSK_Rb250Fd250, +13dbM
  if (!rf95.setFrequency(RF95_FREQ)) {
    Serial.println("setFrequency failed");
    error(3);
    while(1);
  }
  Serial.print("Set Freq to: "); Serial.println(RF95_FREQ);

  // Defaults after init are 434.0MHz, 13dBm, Bw = 125 kHz, Cr = 4/5, Sf = 128chips/symbol, CRC on

  // The default transmitter power is 13dBm, using PA_BOOST.
  // If you are using RFM95/96/97/98 modules which uses the PA_BOOST transmitter pin, then
  // you can set transmitter powers from 5 to 23 dBm:
  rf95.setTxPower(23, false);
  rf95.setPreambleLength(8);
  rf95.setCodingRate4(5);
  // rf95.enableCrc();
  rf95.setSpreadingFactor(8);
  rf95.setSignalBandwidth(125E3);
  // rf95.setSyncWord(SYNC_WORD);
}

int16_t counter = 0;  // packet counter, we increment per xmission
uint8_t buf[RH_RF95_MAX_MESSAGE_LEN];
unsigned long start = 0;

// Device will record values every REPORT_INT and
// attempt to make 30 reports at random intervals
// in that time if it doesn't hear a response
void loop() {
  delay(100);
  unsigned long end = millis();
  unsigned long duration = end - start;

  if (duration / 1000 >= REPORT_INT) {
    start = end;

    // Gymnastics required to get analogRead val into float
    char radiopacket[20] = {' '};
    sprintf(radiopacket, "%s", "hello world");
    
    itoa(counter++, radiopacket+13, 10);
    bool reply_flag = false;

    int attempts = 1;
    while ((!reply_flag) && (attempts < 30) && (((millis() - start) / 1000) < REPORT_INT)) {
      int rand_time = rand() % 5 + 1;
      Serial.print("Attempt "); Serial.println(attempts);
      Serial.print("Sending "); Serial.println(radiopacket);
      digitalWrite(led, HIGH);    // turn the LED on

      delay(10);
    // =============================== BASIC TRANSMISSION ===============================  
      rf95.send((uint8_t *)radiopacket, strlen(radiopacket)+1);

      Serial.println("Waiting for packet to complete...");
      delay(10);
      rf95.waitPacketSent();
      // Now wait for a reply
      uint8_t buf[RH_RF95_MAX_MESSAGE_LEN];
      uint8_t len = sizeof(buf);
      
      Serial.println("Waiting for reply...");
  
      if (rf95.waitAvailableTimeout(5000)) {
        // Should be a reply message for us now
        if (rf95.recv(buf, &len)) {
          Serial.print("Got reply: ");
          Serial.println((char*)buf);
          Serial.print("RSSI: ");
          Serial.println(rf95.lastRssi(), DEC);
          reply_flag = true;
        } else {
          Serial.println("Receive failed");
          break;
        }
      } else {
        Serial.println("No reply, is there a listener around?");
        attempts += 1;
        digitalWrite(led, LOW);    // turn the LED off
        delay(rand_time * 1000);
      }
    // =============================== END BASIC TRANSMISSION ===============================
    }
    digitalWrite(led, LOW);    // turn the LED off
  }
}

// blink out an error code
void error(uint8_t errnum) {
  while(1) {
    uint8_t i;
    for (i=0; i<errnum; i++) {
      digitalWrite(LED_RED, HIGH);
      delay(200);
      digitalWrite(LED_RED, LOW);
      delay(200);
      yield();
    }
    for (i=errnum; i<10; i++) {
      delay(300);
      yield();
    }
  }
}