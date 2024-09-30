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

// Should these be switched to #define?
int led = LED_BUILTIN;
int sense1 = A0; // sensor signal pin 1
int sense2 = A2; // sensor signal pin 2
float val1 = 0.0;  // variable to store the value read
float val2 = 0.0;  // variable to store the value read
int REPORT_INT = 30.0; // interval at which to record and report data

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
  while(!Serial);

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

  // SD Card Initialization
  // see if the card is present and can be initialized:
  if (!SD.begin(SD_CS)) {
    Serial.println("Card init. failed!");
    error(4);
    while(1);
  }
  Serial.println("SD card OK");

  File root = SD.open("/");
  printDirectory(root, 0);
  
  char filename[15];
  strcpy(filename, "/ANALOG00.TXT");
  for (uint8_t i = 0; i < 100; i++) {
    filename[7] = '0' + i/10;
    filename[8] = '0' + i%10;
    // create if does not exist, do not open existing, write, sync after write
    if (! SD.exists(filename)) {
      break;
    }
  }

  logfile = SD.open(filename, FILE_WRITE);
  if( ! logfile ) {
    Serial.print("Couldnt create "); 
    Serial.println(filename);
    error(5);
    while(1);
  }
  Serial.print("Writing to "); 
  Serial.println(filename);
  Serial.println("Ready!");
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
    val1 = analogRead(sense1);
    // Read the data on pin A0 and log it to SD card and to serial
    logfile.print("Sense1 = "); logfile.println(val1);
    Serial.print("Sense1 = "); Serial.println(val1);

    // save the output!
    logfile.flush();

    // Gymnastics required to get analogRead val into float
    char radiopacket[20] = {' '};
    int val1_int = (int) val1;
    float val1_float = (abs(val1) - abs(val1_int)) * 100000;
    int val1_fra = (int)val1_float;
    sprintf(radiopacket, "%d.%d", val1_int, val1_fra);
    
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

    // =============================== DATAGRAM TRANSMISSION ===============================
    // Serial.println("Sending to rf95_reliable_datagram_server");

    // // Send a message to manager_server
    // snprintf((char *)buf, sizeof(buf), "to server counter=%d", ++counter);

    // if (manager.sendtoWait(buf, strlen((char *)buf), SERVER_ADDRESS))
    // {
    // 	// Now wait for a reply from the server
    // 	uint8_t len = sizeof(buf);
    // 	uint8_t from;
    // 	if (manager.recvfromAckTimeout(buf, &len, 2000, &from))
    // 	{
    // 		buf[len] = 0;
    // 		Serial.printf("got reply from 0x%02x rssi=%d %s", from, rf95.lastRssi(), (char *) buf);
    // 	}
    // 	else
    // 	{
    // 		Serial.println("No reply, is rf95_reliable_datagram_server running?");
    // 	}
    // }
    // else
    // 	Serial.println("sendtoWait failed");
    //   error(6);
    // delay(500);
    // =============================== END DATAGRAM TRANSMISSION ===============================
    }
    digitalWrite(led, LOW);    // turn the LED off
  }
}

void printDirectory(File dir, int numTabs) {
   while(true) {
     
     File entry =  dir.openNextFile();
     if (! entry) {
       // no more files
       break;
     }
     for (uint8_t i=0; i<numTabs; i++) {
       Serial.print('\t');
     }
     Serial.print(entry.name());
     if (entry.isDirectory()) {
       Serial.println("/");
       printDirectory(entry, numTabs+1);
     } else {
       // files have sizes, directories do not
       Serial.print("\t\t");
       Serial.println(entry.size(), DEC);
     }
     entry.close();
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