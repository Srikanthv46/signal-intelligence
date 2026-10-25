# 05 — SiK radio (Si443x packet handler) and MAVLink frame details

Research notes, 2026-10-21. Facts only. **[V]** = verified in a primary source (quoted/parsed), **[U]** = not verified, needs a check on a capture.

Sources:
- Si4430 datasheet Prelim. Rev 0.4 (full register descriptions): https://mm.digikey.com/Volume0/opasdata/d220001/medias/docus/2945/Si4430%20Rev0.4.pdf. silabs.com returned HTTP 403 to scripted downloads.
- Si4430/31/32-B1 datasheet Rev 1.1: https://www.nicerf.com/pdf/si4430-31-32.pdf. Figures 23 and 24 were read as images. AN440 itself could not be downloaded (403).
- SiK source (BSD-2): https://github.com/ArduPilot/SiK, files `Firmware/radio/{radio_443x.c,crc.c,mavlink.c,packet.c,tdm.c}` and `Firmware/include/Si1000_defs.h`. A local copy is in `src_sik/`.
- MAVLink c_library_v2 (MIT): https://github.com/mavlink/c_library_v2 (`checksum.h`, `mavlink_helpers.h`, `mavlink_types.h`, `minimal/`, `standard/`, `common/`, `ardupilotmega/`). The tables below were machine-extracted from the generated headers fetched on 2026-10-21.
- PX4 custom mode: https://github.com/PX4/PX4-Autopilot/blob/main/src/modules/commander/px4_custom_mode.h
- A field report of SiK decoding with a Pluto: https://discuss.ardupilot.org/t/encapsulating-mavlink-packets-in-sik-payload/121232

---

## 1. SiK on-air frame (non-Golay, the default `ECC=0`)

| # | Field | Length | Value / source | Whitened? | In HW CRC? |
|---|---|---|---|---|---|
| 1 | Preamble | 16 nibbles = 8 bytes | `1010…` "always a series of alternating ones and zeroes, starting with a one" [V, DS §6.2]; SiK `preamble_length=16` nibbles [V] | no | no |
| 2 | Sync word | 2 bytes | Reg 36h reset value = **0x2D**, reg 37h = **0xD4** [V, register table]. Transmitted in descending order (Sync 3 first) [V, reg 33h synclen]. SiK never writes them [V]. The forum capture shows `2d d4` [V, 3rd-party] | no | no |
| 3 | TX header | 2 bytes | hdlen=010 sends "Header 3 and 2", in descending order [V]. SiK: header3 = netid>>8, header2 = netid&0xFF [V, radio_443x.c:993] | **yes** | yes |
| 4 | Packet length | 1 byte | fixpklen=0 means "the packet length is included in the header" [V, reg 33h]. SiK writes reg 33h = HDLEN_2BYTE\|SYNCLEN_2BYTE, so fixpklen=0 and the length byte **is sent** [V] | **yes** | yes |
| 5 | Data | pklen bytes | SiK payload = serial data (MAVLink frames when framing is on) + TDM trailer (2 bytes) | **yes** | yes |
| 6 | CRC | 2 bytes | CRC-16 (IBM) [V] | **yes** (Fig. 23) | — |

- **Bit order:** MSB first. Reg 30h bit 6 `lsbfrst` = "The LSB of the data will be transmitted/received first if this bit is set" [V]. SiK writes reg 30h = ENPACRX(0x80)\|ENPACTX(0x08)\|ENCRC(0x04)\|CRC_16(0x01) = **0x8D** (Si1000_defs.h) [V], so lsbfrst=0: **MSB first**, and crcdonly=0.
- **Golay mode (`ECC=1`):** reg 30h = 0x88 (no HW CRC) and hdlen=0 [V]. The FIFO then carries Golay(24,12)-encoded `{netid0, netid1, len}`, `{crc_lo, crc_hi, len}` and the payload, with the payload rounded up to a multiple of 3 bytes and the encoded length = (rlen+6)×2 [V, radio_transmit_golay]. This uses the software CRC in §2b.
- **TDM trailer** (last 2 bytes of the data field, tdm.c) [V]: `uint16 window:13, command:1, bonus:1, resend:1`. With AES, a `uint16 crc` follows. SDCC/8051 bitfield order is LSB-first within a little-endian uint16 [U, compiler convention].
- **GFSK frequency polarity** ('1' = +Δf?) [U]. Resolve it on a capture by matching 0x2DD4 against the inverted pattern 0xD22B.

## 2. CRC

### 2a. Hardware CRC (non-Golay)
| Item | Value |
|---|---|
| Polynomial select | reg 30h crc[1:0]: 00 CCITT, **01 CRC-16 (IBM)**, 10 IEC-16, 11 Biacheva [V]. SiK selects 01. |
| Polynomial | CRC-16/IBM = x^16+x^15+x^2+1 = **0x8005** (standard meaning of the name) [U: the datasheet names it but does not give the hex] |
| Coverage | crcdonly (reg 30h bit 5) = "When this bit is set to 1 the CRC is calculated on and checked against the packet data fields only" [V]. SiK leaves it 0, so **CRC covers header(2) + length(1) + data** (Rev 1.1 §6.5: "CRC can be applied to only the data portion of the packet or to the data, packet length and header fields") [V]. Note: this is *not* data-only. |
| Whitening | The CRC field lies inside the whitening span in Fig. 23 [V], so the CRC is computed on **unwhitened** bytes and then whitened along with them (the usual reading of the figure) [U on hardware] |
| Seed / init | **Not stated** in either datasheet revision [U]. Reg 30h has no seed bit (Rev 1.1 bit 4 = `skip2ph`, Rev 0.4 = reserved). The later Si446x has a CRC_SEED bit (0 = all 0s, 1 = all 1s) (AN626), which suggests 0x0000 or 0xFFFF. |
| Reflection, xorout, byte order | **Not stated** [U]. Given MSB-first transmission, the likely form is non-reflected with the CRC high byte first, but this is unverified. |
| Practical test | On the first decoded frames, brute-force {init 0x0000, 0xFFFF} × {refin/refout} × {xorout 0, 0xFFFF} × {byte order}, over header+len+data. The forum (3rd party) reports that "standard CRCs" did not match, unresolved, so expect to have to search. |

### 2b. SiK software crc16 (Golay mode, AES trailer.crc, fhop seed) [V, crc.c]
- Table-driven, `high = low = 0`, then each byte is shifted into `low`.
- I ran the actual tables from crc.c in Python: the result equals the **non-augmented remainder M(x) mod 0x1021 with init 0**. This is not CRC-16/XMODEM (no 16-bit augmentation).
- crc16("123456789") = **0xBEEF**.
- On air: byte0 = crc&0xFF, byte1 = crc>>8.

## 3. Data whitening (reg 70h bit 0 `enwhite`)
| Item | Value |
|---|---|
| Generator | "payload data bits are XORed with a pseudorandom sequence output from the built-in PN9 generator. The generator is initialized at the beginning of the payload" [V, §6.5] |
| Span | Fig. 23 (both revisions): whitening covers **Header/Address + PK Length + Data + CRC**, and not the preamble or sync [V]. "Payload" therefore means everything after the sync word. Reset per packet [V, "initialized at the beginning"]. |
| Polynomial | PN9 = x^9+x^5+1 (generic PN9 definition) [U for Si443x specifically] |
| Seed | 0x1FF (all ones) [U]. A 3rd-party Pluto/URH user reported dewhitening "starting at the Net ID field" with PN9 seed 0x1FF recovering the netid/length/MAVLink header. |
| Bit/byte order | [U]. Starting hypothesis: the CC1101-style PN9 byte sequence from 0x1FF (FF E1 1D 9A ED 85 33 24 …) XORed MSB-first per byte. Verify: the dewhitened header must equal the netid and the length must match the burst duration. |
| Packet-handler-off note | "Manchester encoding is supported but data whitening, CRC, or header checks [are] not" when the packet handler is disabled [V] |

## 4. Manchester (reg 70h) [V]
- `enmanch` (bit 1): Manchester coding on. The chip rate doubles (preamble etc.) and the effective rate is limited to 128 kbps (Rev 1.1) or 64 kbps (Rev 0.4).
- Manchester spans the whole packet (preamble through CRC) [V, Fig. 23].
- `enmaninv` (bit 2): "Manchester Data Inversion is Enabled". Fig. 24 shows it swaps the two half-bit symbols, i.e. inverts the encoded chips.
- `manppol` (bit 3): the TX preamble is "a series of 1 if set, or series of 0 if reset". It affects only the TX side and is valid only with Manchester on.
- SiK default 0x0D/0x2D means enmanch=0, so enmaninv and manppol should have no effect [U: there is no datasheet statement on enmaninv with enmanch=0, so check the polarity on a capture].
- Reg 70h reset = 0x0C [V]. Bit 5 `txdtrtscale` is "set for Data Rates below 30 kbps" (hence 0x2D) [V].

## 5. MAVLink framing [V, c_library_v2]

**v1:** `FE | len | seq | sysid | compid | msgid(1) | payload[len] | crc_lo | crc_hi`
- Header 6 bytes, total 8+len.
- v1 sends the full base payload (`min_length`, no extensions).

**v2:** `FD | len | incompat_flags | compat_flags | seq | sysid | compid | msgid(3, LE) | payload[len] | crc_lo | crc_hi | [signature 13]`
- Header 10 bytes, total 12+len (+13 if `incompat_flags & 0x01`, MAVLINK_IFLAG_SIGNED).
- **Trailing zero bytes of the payload are trimmed** (`_mav_trim_payload`). The receiver zero-fills up to the full length, so `len` can be less than the minimum length.

**Checksum:** CRC-16/MCRF4XX.
- `X25_INIT_CRC 0xffff`, byte-wise `tmp=b^(crc&0xff); tmp^=tmp<<4; crc=(crc>>8)^(tmp<<8)^(tmp<<3)^(tmp>>4)`. This is poly 0x1021 reflected (0x8408), no xorout.
- Computed over **header bytes after STX through the end of the payload, then one extra byte CRC_EXTRA** (`crc_calculate(&buf[1], header_len-1)` + payload + `crc_accumulate(crc_extra)`).
- Sent little-endian.
- All multi-byte fields are little-endian. Fields are sorted by type size (8 > 4 > 2 > 1 byte, stable within a size). Extension fields are appended after the base fields, in XML order.

### Message table
Copied from `MAVLINK_MSG_ID_xxx_{MIN_LEN,LEN,CRC}` and the MESSAGE_INFO offsets.
- Format of the "Wire layout" column: `name type@byte_offset`.
- Fields after **ext:** are MAVLink-2 extensions (absent in v1).
- Header directories: minimal/ (0), standard/ (33, 148), ardupilotmega/ (166), common/ (all others).

| ID | Name | min/max len | CRC_EXTRA | Wire layout |
|---|---|---|---|---|
| 0 | HEARTBEAT | 9 / 9 | 50 | custom_mode u32@0; type u8@4; autopilot u8@5; base_mode u8@6; system_status u8@7; mavlink_version u8@8 |
| 1 | SYS_STATUS | 31 / 43 | 124 | onboard_control_sensors_present u32@0; _enabled u32@4; _health u32@8; load u16@12; voltage_battery u16@14; current_battery i16@16; drop_rate_comm u16@18; errors_comm u16@20; errors_count1..4 u16@22,24,26,28; battery_remaining i8@30; **ext:** sensors_present_extended u32@31; _enabled_extended u32@35; _health_extended u32@39 |
| 2 | SYSTEM_TIME | 12 / 12 | 137 | time_unix_usec u64@0; time_boot_ms u32@8 |
| 22 | PARAM_VALUE | 25 / 25 | 220 | param_value f32@0; param_count u16@4; param_index u16@6; param_id char[16]@8; param_type u8@24 |
| 24 | GPS_RAW_INT | 30 / 52 | 24 | time_usec u64@0; lat i32@8; lon i32@12; alt i32@16; eph u16@20; epv u16@22; vel u16@24; cog u16@26; fix_type u8@28; satellites_visible u8@29; **ext:** alt_ellipsoid i32@30; h_acc u32@34; v_acc u32@38; vel_acc u32@42; hdg_acc u32@46; yaw u16@50 |
| 27 | RAW_IMU | 26 / 29 | 144 | time_usec u64@0; xacc,yacc,zacc i16@8,10,12; xgyro,ygyro,zgyro i16@14,16,18; xmag,ymag,zmag i16@20,22,24; **ext:** id u8@26; temperature i16@27 |
| 29 | SCALED_PRESSURE | 14 / 16 | 115 | time_boot_ms u32@0; press_abs f32@4; press_diff f32@8; temperature i16@12; **ext:** temperature_press_diff i16@14 |
| 30 | ATTITUDE | 28 / 28 | 39 | time_boot_ms u32@0; roll f32@4; pitch f32@8; yaw f32@12; rollspeed f32@16; pitchspeed f32@20; yawspeed f32@24 |
| 33 | GLOBAL_POSITION_INT | 28 / 28 | 104 | time_boot_ms u32@0; lat i32@4; lon i32@8; alt i32@12; relative_alt i32@16; vx i16@20; vy i16@22; vz i16@24; hdg u16@26 |
| 35 | RC_CHANNELS_RAW | 22 / 22 | 244 | time_boot_ms u32@0; chan1..8_raw u16@4..18; port u8@20; rssi u8@21 |
| 36 | SERVO_OUTPUT_RAW | 21 / 37 | 222 | time_usec u32@0; servo1..8_raw u16@4..18; port u8@20; **ext:** servo9..16_raw u16@21..35 |
| 42 | MISSION_CURRENT | 2 / 18 | 28 | seq u16@0; **ext:** total u16@2; mission_state u8@4; mission_mode u8@5; mission_id u32@6; fence_id u32@10; rally_points_id u32@14 |
| 62 | NAV_CONTROLLER_OUTPUT | 26 / 26 | 183 | nav_roll f32@0; nav_pitch f32@4; alt_error f32@8; aspd_error f32@12; xtrack_error f32@16; nav_bearing i16@20; target_bearing i16@22; wp_dist u16@24 |
| 65 | RC_CHANNELS | 42 / 42 | 118 | time_boot_ms u32@0; chan1..18_raw u16@4..38; chancount u8@40; rssi u8@41 |
| 74 | VFR_HUD | 20 / 20 | 20 | airspeed f32@0; groundspeed f32@4; alt f32@8; climb f32@12; heading i16@16; throttle u16@18 |
| 76 | COMMAND_LONG | 33 / 33 | 152 | param1..7 f32@0..24; command u16@28; target_system u8@30; target_component u8@31; confirmation u8@32 |
| 77 | COMMAND_ACK | 3 / 10 | 143 | command u16@0; result u8@2; **ext:** progress u8@3; result_param2 i32@4; target_system u8@8; target_component u8@9 |
| 109 | RADIO_STATUS | 9 / 9 | 185 | rxerrors u16@0; fixed u16@2; rssi u8@4; remrssi u8@5; txbuf u8@6; noise u8@7; remnoise u8@8 |
| 111 | TIMESYNC | 16 / 18 | 34 | tc1 i64@0; ts1 i64@8; **ext:** target_system u8@16; target_component u8@17 |
| 147 | BATTERY_STATUS | 36 / 54 | 154 | current_consumed i32@0; energy_consumed i32@4; temperature i16@8; voltages u16[10]@10; current_battery i16@30; id u8@32; battery_function u8@33; type u8@34; battery_remaining i8@35; **ext:** time_remaining i32@36; charge_state u8@40; voltages_ext u16[4]@41; mode u8@49; fault_bitmask u32@50 |
| 148 | AUTOPILOT_VERSION | 60 / 78 | 178 | capabilities u64@0; uid u64@8; flight_sw_version u32@16; middleware_sw_version u32@20; os_sw_version u32@24; board_version u32@28; vendor_id u16@32; product_id u16@34; flight_custom_version u8[8]@36; middleware_custom_version u8[8]@44; os_custom_version u8[8]@52; **ext:** uid2 u8[18]@60 |
| 166 | RADIO (ardupilotmega) | 9 / 9 | 21 | same layout as RADIO_STATUS: rxerrors u16@0; fixed u16@2; rssi u8@4; remrssi u8@5; txbuf u8@6; noise u8@7; remnoise u8@8 |
| 242 | HOME_POSITION | 52 / 60 | 104 | latitude i32@0; longitude i32@4; altitude i32@8; x f32@12; y f32@16; z f32@20; q f32[4]@24; approach_x f32@40; approach_y f32@44; approach_z f32@48; **ext:** time_usec u64@52 |
| 245 | EXTENDED_SYS_STATE | 2 / 2 | 130 | vtol_state u8@0; landed_state u8@1 |
| 253 | STATUSTEXT | 51 / 54 | 83 | severity u8@0; text char[50]@1; **ext:** id u16@51; chunk_seq u8@53 |

Units, per the XML: lat/lon degE7; alt mm; relative_alt mm; GLOBAL_POSITION_INT vx/vy/vz cm/s; hdg cdeg (65535 = unknown); ATTITUDE rad and rad/s; VFR_HUD heading deg, throttle %; GPS_RAW_INT eph/epv ×100, vel cm/s, cog cdeg.

### SiK-generated MAVLink [V, SiK `mavlink.c`, `packet.c`]
- `MAVLink_report()` builds a hand-crafted **RADIO_STATUS (109, CRC_EXTRA 185)** with sysid `'3'` (0x33) and compid `'D'` (0x44).
  - v1 or v2 framing, matching what it has seen from the autopilot. In v2: incompat=compat=0 and msgid bytes = 109,0,0.
  - Payload: rxerrors=errors.rx_errors, fixed=corrected_packets, txbuf=serial_read_space(), rssi/remrssi=average RSSI (local/remote), noise/remnoise.
- On a PROTOCOL_VERSION handshake it sends **PROTOCOL_VERSION (300, CRC_EXTRA 217)** instead.
- **It is written to the local serial port (`serial_write_buf`), not to the air.** So sysid 0x33 / compid 0x44 RADIO_STATUS frames are not expected on air [V by code reading].
- `seen_mavlink` detection: a v1 heartbeat is `buf[0]==0xFE && buf[1]==9 && buf[5]==0`; a v2 heartbeat is `0xFD, len<=9, msgid bytes 0`.
- With MAVLink framing (`MAVLINK=1`, the default) `mavlink_frame()` packs **whole MAVLink frames** into each radio packet when they fit. On air, the data field therefore usually starts with 0xFE/0xFD [V, packet.c].

## 6. HEARTBEAT enums [V, minimal.h / ardupilotmega.h]

**MAV_TYPE:**
0 GENERIC, 1 FIXED_WING, 2 QUADROTOR, 3 COAXIAL, 4 HELICOPTER, 5 ANTENNA_TRACKER, 6 GCS, 7 AIRSHIP, 8 FREE_BALLOON, 9 ROCKET, 10 GROUND_ROVER, 11 SURFACE_BOAT, 12 SUBMARINE, 13 HEXAROTOR, 14 OCTOROTOR, 15 TRICOPTER, 16 FLAPPING_WING, 17 KITE, 18 ONBOARD_CONTROLLER, 19 VTOL_TAILSITTER_DUOROTOR, 20 VTOL_TAILSITTER_QUADROTOR, 21 VTOL_TILTROTOR, 22 VTOL_FIXEDROTOR, 23 VTOL_TAILSITTER, 24 VTOL_TILTWING, 25 VTOL_RESERVED5, 26 GIMBAL, 27 ADSB, 28 PARAFOIL, 29 DODECAROTOR, 30 CAMERA, 31 CHARGING_STATION, 32 FLARM, 33 SERVO, 34 ODID, 35 DECAROTOR, 36 BATTERY, 37 PARACHUTE, 38 LOG, 39 OSD, 40 IMU, 41 GPS, 42 WINCH, 43 GENERIC_MULTIROTOR, 44 ILLUMINATOR, 45 SPACECRAFT_ORBITER, 46 GROUND_QUADRUPED, 47 VTOL_GYRODYNE, 48 GRIPPER, 49 RADIO.

**MAV_AUTOPILOT:**
0 GENERIC, 1 RESERVED, 2 SLUGS, **3 ARDUPILOTMEGA**, 4 OPENPILOT, 5 GENERIC_WAYPOINTS_ONLY, 6 GENERIC_WAYPOINTS_AND_SIMPLE_NAVIGATION_ONLY, 7 GENERIC_MISSION_FULL, 8 INVALID (used by GCS/peripherals), 9 PPZ, 10 UDB, 11 FP, **12 PX4**, 13 SMACCMPILOT, 14 AUTOQUAD, 15 ARMAZILA, 16 AEROB, 17 ASLUAV, 18 SMARTAP, 19 AIRRAILS, 20 REFLEX, 21 FLIX.

**base_mode (MAV_MODE_FLAG bits):**

| Bit value | Flag |
|---|---|
| 0x01 | CUSTOM_MODE_ENABLED |
| 0x02 | TEST_ENABLED |
| 0x04 | AUTO_ENABLED |
| 0x08 | GUIDED_ENABLED |
| 0x10 | STABILIZE_ENABLED |
| 0x20 | HIL_ENABLED |
| 0x40 | MANUAL_INPUT_ENABLED |
| 0x80 | SAFETY_ARMED |

**system_status (MAV_STATE):** 0 UNINIT, 1 BOOT, 2 CALIBRATING, 3 STANDBY, 4 ACTIVE, 5 CRITICAL, 6 EMERGENCY, 7 POWEROFF, 8 FLIGHT_TERMINATION.

**ArduPilot custom_mode** (autopilot=3). Pick the table by MAV_TYPE.

| Vehicle | custom_mode values |
|---|---|
| Copter (COPTER_MODE) | 0 STABILIZE, 1 ACRO, 2 ALT_HOLD, 3 AUTO, 4 GUIDED, 5 LOITER, 6 RTL, 7 CIRCLE, 9 LAND, 11 DRIFT, 13 SPORT, 14 FLIP, 15 AUTOTUNE, 16 POSHOLD, 17 BRAKE, 18 THROW, 19 AVOID_ADSB, 20 GUIDED_NOGPS, 21 SMART_RTL, 22 FLOWHOLD, 23 FOLLOW, 24 ZIGZAG, 25 SYSTEMID, 26 AUTOROTATE, 27 AUTO_RTL, 28 TURTLE |
| Plane (PLANE_MODE) | 0 MANUAL, 1 CIRCLE, 2 STABILIZE, 3 TRAINING, 4 ACRO, 5 FBWA, 6 FBWB, 7 CRUISE, 8 AUTOTUNE, 10 AUTO, 11 RTL, 12 LOITER, 13 TAKEOFF, 14 AVOID_ADSB, 15 GUIDED, 16 INITIALIZING, 17 QSTABILIZE, 18 QHOVER, 19 QLOITER, 20 QLAND, 21 QRTL, 22 QAUTOTUNE, 23 QACRO, 24 THERMAL, 25 LOITER_ALT_QLAND, 26 AUTOLAND |
| Rover (ROVER_MODE) | 0 MANUAL, 1 ACRO, 3 STEERING, 4 HOLD, 5 LOITER, 6 FOLLOW, 7 SIMPLE, 8 DOCK, 9 CIRCLE, 10 AUTO, 11 RTL, 12 SMART_RTL, 15 GUIDED, 16 INITIALIZING |

**PX4 custom_mode** (autopilot=12) [V, px4_custom_mode.h]: a union of `{uint16 reserved; uint8 main_mode; uint8 sub_mode}`. As a little-endian u32: main_mode = (custom_mode>>16)&0xFF, sub_mode = custom_mode>>24.
- main_mode: 1 MANUAL, 2 ALTCTL, 3 POSCTL, 4 AUTO, 5 ACRO, 6 OFFBOARD, 7 STABILIZED, 8 RATTITUDE_LEGACY, 9 SIMPLE, 10 TERMINATION, 11 ALTITUDE_CRUISE, 12 MANUAL_PARKING.
- AUTO sub_mode: 1 READY, 2 TAKEOFF, 3 LOITER, 4 MISSION, 5 RTL, 6 LAND, 7 reserved (was RTGS), 8 FOLLOW_TARGET, 9 PRECLAND, 10 VTOL_TAKEOFF, 11–18 EXTERNAL1–8, 19 GUIDED_COURSE, 20 AUTO_DESCEND.
- POSCTL sub_mode: 0 POSCTL, 1 ORBIT, 2 SLOW.

## 7. Open items to settle on a capture
1. HW CRC-16 seed, reflection and byte order (§2a).
2. The exact PN9 byte sequence and bit order for Si443x whitening (§3).
3. GFSK polarity, and whether enmaninv=1 has any effect with enmanch=0.
4. TDM trailer bit layout.

A known-content capture solves all four:
- a SiK pair at a known NETID;
- a known MAVLink stream, e.g. HEARTBEAT with known sysid, where CRC_EXTRA lets the MAVLink layer be checked independently of the radio CRC.
