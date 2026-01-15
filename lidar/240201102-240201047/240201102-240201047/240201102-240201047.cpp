#define _CRT_SECURE_NO_WARNINGS
#define NOMINMAX

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <random>
#include <limits>
#include <algorithm>
#include <utility>
#include <tuple>
#include <iomanip>

#include <Windows.h>
#include <d3d11.h>
#include <tchar.h>

#pragma comment(lib, "d3d11.lib")
#pragma comment(lib, "dxgi.lib")

#include "lib/curl-msvc/include/curl/curl.h"
#pragma comment(linker, "/LIBPATH:\"lib\\curl-msvc\\lib\\dll-release-x64\"")
#pragma comment(lib, "libcurl.dll")
#pragma comment(lib, "ws2_32.lib")
#pragma comment(lib, "crypt32.lib")
#pragma comment(lib, "wldap32.lib")

#include "lib/imgui.h"
#include "lib/imgui_impl_win32.h"
#include "lib/imgui_impl_dx11.h"
#include "lib/implot.h"

using namespace std;

struct Koordinat {
	double k_x = 0.0;
	double k_y = 0.0;
};

struct DogruModeli {
	double kat_A = 0.0;
	double kat_B = -1.0;
	double sabit_C = 0.0;
};

struct SegmentBilgisi {
	DogruModeli model;
	vector<Koordinat> ait_noktalar;
	Koordinat uc_1;
	Koordinat uc_2;
	double mesafe = 0.0;
	bool gorunur = true;
};

size_t VeriyiDosyayaYaz(void* gelen_veri_blogu, size_t eleman_boyutu, size_t eleman_sayisi, FILE* hedef_dosya) {
	size_t yazilan_boyut = fwrite(gelen_veri_blogu, eleman_boyutu, eleman_sayisi, hedef_dosya);
	return yazilan_boyut;
}

bool YapilandirmaDosyasiniIndir(const char* kaynak_url, const char* yerel_dosya_adi) {
	CURL* curl_oturum;
	FILE* dosya_akisi;
	CURLcode sonuc_kodu;

	curl_global_init(CURL_GLOBAL_DEFAULT);
	curl_oturum = curl_easy_init();

	if (curl_oturum) {
		dosya_akisi = fopen(yerel_dosya_adi, "wb");
		if (dosya_akisi == NULL) {
			MessageBoxA(NULL, "Yerel config.toml dosyasi olusturulamadi/acilamadi.", "Dosya Yazma Hatasi", MB_OK | MB_ICONERROR);
			curl_easy_cleanup(curl_oturum);
			curl_global_cleanup();
			return false;
		}

		curl_easy_setopt(curl_oturum, CURLOPT_URL, kaynak_url);
		curl_easy_setopt(curl_oturum, CURLOPT_WRITEFUNCTION, VeriyiDosyayaYaz);
		curl_easy_setopt(curl_oturum, CURLOPT_WRITEDATA, dosya_akisi);

		sonuc_kodu = curl_easy_perform(curl_oturum);

		fclose(dosya_akisi);
		curl_easy_cleanup(curl_oturum);

		if (sonuc_kodu != CURLE_OK) {
			string hata_metni = "config.toml indirilemedi. libcurl hatasi: ";
			hata_metni += curl_easy_strerror(sonuc_kodu);
			MessageBoxA(NULL, hata_metni.c_str(), "Indirme Hatasi", MB_OK | MB_ICONERROR);
			curl_global_cleanup();
			return false;
		}
	}
	else {
		MessageBoxA(NULL, "libcurl baslatilamadi (curl_easy_init failed).", "libcurl Hatasi", MB_OK | MB_ICONERROR);
		curl_global_cleanup();
		return false;
	}

	curl_global_cleanup();
	return true;
}

vector<Koordinat> KutupsalDonusum(const vector<float>& olcumler,
	float baslangic_radyan, float artis_radyan, float min_mesafe, float max_mesafe) {

	vector<Koordinat> koordinatlar;

	for (size_t i = 0; i < olcumler.size(); ++i) {
		float mesafe = olcumler[i];

		if (mesafe == -1.0f) continue;
		if (mesafe < min_mesafe || mesafe > max_mesafe) continue;
		if (isinf(mesafe) || isnan(mesafe)) continue;

		float radyan = baslangic_radyan + i * artis_radyan;
		float gecici_x = mesafe * cos(radyan);
		float gecici_y = mesafe * sin(radyan);

		Koordinat nokta;
		nokta.k_x = gecici_x;
		nokta.k_y = gecici_y;
		koordinatlar.push_back(nokta);
	}
	return koordinatlar;
}

DogruModeli NoktalardanModelOlustur(const Koordinat& nokta_A, const Koordinat& nokta_B) {
	DogruModeli yeni_model;

	if (fabs(nokta_B.k_x - nokta_A.k_x) < 1e-6) {
		yeni_model.kat_A = numeric_limits<double>::infinity();
		yeni_model.sabit_C = nokta_A.k_x;
	}
	else {
		yeni_model.kat_A = (nokta_B.k_y - nokta_A.k_y) / (nokta_B.k_x - nokta_A.k_x);
		yeni_model.sabit_C = nokta_A.k_y - yeni_model.kat_A * nokta_A.k_x;
	}
	return yeni_model;
}

double NoktaninModeleUzakligi(const Koordinat& nokta, const DogruModeli& model) {
	if (isinf(model.kat_A)) {
		return fabs(nokta.k_x - model.sabit_C);
	}

	return fabs(model.kat_A * nokta.k_x + model.kat_B * nokta.k_y + model.sabit_C) / sqrt(model.kat_A * model.kat_A + 1);
}

pair<Koordinat, Koordinat> SegmentUclariBul(const vector<Koordinat>& nokta_listesi, const DogruModeli& model) {
	if (nokta_listesi.size() < 2) {
		return { Koordinat(), Koordinat() };
	}

	Koordinat uc_A = nokta_listesi[0];
	Koordinat uc_B = nokta_listesi[0];

	if (isinf(model.kat_A)) {
		for (const auto& gecerli_nokta : nokta_listesi) {
			if (gecerli_nokta.k_y < uc_A.k_y) uc_A = gecerli_nokta;
			if (gecerli_nokta.k_y > uc_B.k_y) uc_B = gecerli_nokta;
		}
	}
	else {
		for (const auto& gecerli_nokta : nokta_listesi) {
			if (gecerli_nokta.k_x < uc_A.k_x) uc_A = gecerli_nokta;
			if (gecerli_nokta.k_x > uc_B.k_x) uc_B = gecerli_nokta;
		}
	}
	return { uc_A, uc_B };
}

vector<SegmentBilgisi> RansacModelBulucu(const vector<Koordinat>& nokta_bulutu) {
	vector<SegmentBilgisi> bulunan_segmentler;
	if (nokta_bulutu.size() < 8) return bulunan_segmentler;

	mt19937 rastgele_uretici(17);
	uniform_int_distribution<> indeks_dagilimi(0, static_cast<int>(nokta_bulutu.size()) - 1);

	size_t N = nokta_bulutu.size();

	vector<vector<bool>> denenen_ciftler(N, vector<bool>(N, false));
	const double mesafe_toleransi = 0.01;
	const int min_nokta_sayisi = 8;
	const double birlestirme_toleransi = 0.05;

	bool dongu_devam = true;
	while (dongu_devam) {

		int indeks_1 = indeks_dagilimi(rastgele_uretici);
		int indeks_2 = indeks_dagilimi(rastgele_uretici);
		if (indeks_1 == indeks_2) continue;

		int i_kucuk = min(indeks_1, indeks_2);
		int i_buyuk = max(indeks_1, indeks_2);

		if (denenen_ciftler[i_kucuk][i_buyuk]) continue;

		denenen_ciftler[i_kucuk][i_buyuk] = true;

		const Koordinat& nokta_1 = nokta_bulutu[indeks_1];
		const Koordinat& nokta_2 = nokta_bulutu[indeks_2];

		if (sqrt(pow(nokta_1.k_x - nokta_2.k_x, 2) + pow(nokta_1.k_y - nokta_2.k_y, 2)) < 1e-6)
			continue;

		DogruModeli hipotez_model = NoktalardanModelOlustur(nokta_1, nokta_2);

		vector<Koordinat> destekci_noktalar;
		for (const auto& test_noktasi : nokta_bulutu) {
			double uzaklik = NoktaninModeleUzakligi(test_noktasi, hipotez_model);
			if (uzaklik < mesafe_toleransi) {
				destekci_noktalar.push_back(test_noktasi);
			}
		}

		if (destekci_noktalar.size() >= static_cast<size_t>(min_nokta_sayisi)) {

			double ort_x = 0, ort_y = 0;
			for (auto& nokta : destekci_noktalar) {
				ort_x += nokta.k_x;
				ort_y += nokta.k_y;
			}
			ort_x /= destekci_noktalar.size();
			ort_y /= destekci_noktalar.size();

			double pay = 0, payda = 0;
			for (auto& nokta : destekci_noktalar) {
				pay += (nokta.k_x - ort_x) * (nokta.k_y - ort_y);
				payda += (nokta.k_x - ort_x) * (nokta.k_x - ort_x);
			}

			DogruModeli iyilestirilmis_model;
			if (fabs(payda) < 1e-9) {
				iyilestirilmis_model.kat_A = numeric_limits<double>::infinity();
				iyilestirilmis_model.sabit_C = ort_x;
			}
			else {
				double egim = pay / payda;
				double kesisim_c = ort_y - egim * ort_x;
				iyilestirilmis_model.kat_A = egim;
				iyilestirilmis_model.sabit_C = kesisim_c;
			}

			vector<size_t> nihai_indeksler;
			vector<Koordinat> nihai_noktalar;
			for (size_t j = 0; j < N; ++j) { // N, fonksiyonun başındaki 'nokta_bulutu.size()'
				double uzaklik = NoktaninModeleUzakligi(nokta_bulutu[j], iyilestirilmis_model);
				if (uzaklik < mesafe_toleransi) {
					nihai_indeksler.push_back(j);
					nihai_noktalar.push_back(nokta_bulutu[j]);
				}
			}

			for (size_t i = 0; i < nihai_indeksler.size(); ++i) {
				for (size_t j = i + 1; j < nihai_indeksler.size(); ++j) {
					int idx_a = static_cast<int>(nihai_indeksler[i]);
					int idx_b = static_cast<int>(nihai_indeksler[j]);

					denenen_ciftler[min(idx_a, idx_b)][max(idx_a, idx_b)] = true;
				}
			}

			bool benzer_bulundu = false;
			for (auto& kayitli_segment : bulunan_segmentler) {
				const auto& kayitli_model = kayitli_segment.model;

				if (isinf(iyilestirilmis_model.kat_A) && isinf(kayitli_model.kat_A)) {
					if (fabs(iyilestirilmis_model.sabit_C - kayitli_model.sabit_C) < birlestirme_toleransi) {
						benzer_bulundu = true;
						break;
					}
				}
				else if (!isinf(iyilestirilmis_model.kat_A) && !isinf(kayitli_model.kat_A)) {
					if ((fabs(iyilestirilmis_model.kat_A - kayitli_model.kat_A) < 0.05) && (fabs(iyilestirilmis_model.sabit_C - kayitli_model.sabit_C) < birlestirme_toleransi)) {
						benzer_bulundu = true;
						break;
					}
				}
			}
			if (benzer_bulundu) continue;

			auto [segment_bas, segment_son] = SegmentUclariBul(nihai_noktalar, iyilestirilmis_model);

			double segment_uzunlugu = sqrt(pow(segment_bas.k_x - segment_son.k_x, 2) + pow(segment_bas.k_y - segment_son.k_y, 2));

			SegmentBilgisi yeni_segment;
			yeni_segment.model = iyilestirilmis_model;
			yeni_segment.ait_noktalar = nihai_noktalar;
			yeni_segment.uc_1 = segment_bas;
			yeni_segment.uc_2 = segment_son;
			yeni_segment.mesafe = segment_uzunlugu;
			bulunan_segmentler.push_back(yeni_segment);
		}

		bool tum_ciftler_denendi = true;
		for (size_t i = 0; i < N; ++i) {
			for (size_t j = i + 1; j < N; ++j) {
				if (!denenen_ciftler[i][j]) {
					tum_ciftler_denendi = false;
					break;
				}
			}
			if (!tum_ciftler_denendi) {
				break;
			}
		}

		if (tum_ciftler_denendi) dongu_devam = false;

	}

	return bulunan_segmentler;
}

vector<SegmentBilgisi> SegmentleriBoslugaGoreBol(const vector<Koordinat>& tum_nokta_bulutu) {

	const int min_parca_nokta_sayisi = 8;
	const double max_nokta_arasi_mesafe = 0.50;
	vector<SegmentBilgisi> kaba_segmentler = RansacModelBulucu(tum_nokta_bulutu);
	vector<SegmentBilgisi> ayrilmis_segmentler;

	for (const auto& kaba_segment : kaba_segmentler) {

		vector<Koordinat> siralanacak_noktalar = kaba_segment.ait_noktalar;

		if (siralanacak_noktalar.size() < static_cast<size_t>(min_parca_nokta_sayisi)) {
			continue;
		}

		const double model_egimi = kaba_segment.model.kat_A;

		sort(siralanacak_noktalar.begin(), siralanacak_noktalar.end(), [model_egimi](const Koordinat& a, const Koordinat& b) {
			if (isinf(model_egimi) || fabs(model_egimi) > 1.0) {
				return a.k_y < b.k_y;
			}
			else {
				return a.k_x < b.k_x;
			}
			});

		vector<Koordinat> gecici_parca_noktalari;
		gecici_parca_noktalari.push_back(siralanacak_noktalar[0]);

		for (size_t i = 0; i < siralanacak_noktalar.size() - 1; ++i) {
			const Koordinat& nokta_onceki = siralanacak_noktalar[i];
			const Koordinat& nokta_sonraki = siralanacak_noktalar[i + 1];

			double aradaki_bosluk = sqrt(pow(nokta_onceki.k_x - nokta_sonraki.k_x, 2) + pow(nokta_onceki.k_y - nokta_sonraki.k_y, 2));

			if (aradaki_bosluk > max_nokta_arasi_mesafe) {

				if (gecici_parca_noktalari.size() >= static_cast<size_t>(min_parca_nokta_sayisi)) {
					SegmentBilgisi yeni_parca;
					yeni_parca.model = kaba_segment.model;
					yeni_parca.ait_noktalar = gecici_parca_noktalari;

					auto [parca_basi, parca_sonu] = SegmentUclariBul(gecici_parca_noktalari, yeni_parca.model);
					yeni_parca.uc_1 = parca_basi;
					yeni_parca.uc_2 = parca_sonu;
					yeni_parca.mesafe = sqrt(pow(parca_basi.k_x - parca_sonu.k_x, 2) + pow(parca_basi.k_y - parca_sonu.k_y, 2));

					ayrilmis_segmentler.push_back(yeni_parca);
				}

				gecici_parca_noktalari.clear();
				gecici_parca_noktalari.push_back(nokta_sonraki);
			}
			else {
				gecici_parca_noktalari.push_back(nokta_sonraki);
			}
		}

		if (gecici_parca_noktalari.size() >= static_cast<size_t>(min_parca_nokta_sayisi)) {
			SegmentBilgisi son_parca;
			son_parca.model = kaba_segment.model;
			son_parca.ait_noktalar = gecici_parca_noktalari;

			auto [parca_basi, parca_sonu] = SegmentUclariBul(gecici_parca_noktalari, son_parca.model);
			son_parca.uc_1 = parca_basi;
			son_parca.uc_2 = parca_sonu;
			son_parca.mesafe = sqrt(pow(parca_basi.k_x - parca_sonu.k_x, 2) + pow(parca_basi.k_y - parca_sonu.k_y, 2));

			ayrilmis_segmentler.push_back(son_parca);
		}
	}

	return ayrilmis_segmentler;
}

bool ModelKesisimBul(const DogruModeli& model_1, const DogruModeli& model_2, Koordinat& kesisim_noktasi) {

	bool m1_dikey = isinf(model_1.kat_A);
	bool m2_dikey = isinf(model_2.kat_A);

	// 1. Durum: İkisi de dikeyse
	if (m1_dikey && m2_dikey) {
		return false; // İki dikey doğru her zaman paraleldir.
	}
	// 2. Durum: Sadece model_1 dikeyse (x = k)
	else if (m1_dikey) {
		// x = model_1.sabit_C
		kesisim_noktasi.k_x = model_1.sabit_C;
		// A2*x + B2*y + C2 = 0
		// A2*k + B2*y + C2 = 0
		// B2*y = -C2 - A2*k
		// y = (-C2 - A2*k) / B2
		// Bizim B2 hep -1 olduğu için: y = C2 + A2*k
		kesisim_noktasi.k_y = model_2.sabit_C + model_2.kat_A * model_1.sabit_C;
		return true;
	}
	// 3. Durum: Sadece model_2 dikeyse (x = k)
	else if (m2_dikey) {
		kesisim_noktasi.k_x = model_2.sabit_C;
		// Benzer şekilde: y = C1 + A1*k
		kesisim_noktasi.k_y = model_1.sabit_C + model_1.kat_A * model_2.sabit_C;
		return true;
	}

	// 4. Durum: İkisi de dikey değilse (Orijinal Cramer Kuralı)
	double determinant = model_1.kat_A * model_2.kat_B - model_1.kat_B * model_2.kat_A;

	// Paralel kontrolü
	if (fabs(determinant) < 1e-9) return false;

	kesisim_noktasi.k_x = (model_2.kat_B * (-model_1.sabit_C) - model_1.kat_B * (-model_2.sabit_C)) / determinant;
	kesisim_noktasi.k_y = (model_1.kat_A * (-model_2.sabit_C) - model_2.kat_A * (-model_1.sabit_C)) / determinant;
	return true;
}

bool NoktaSegmentUzerindeMi(const Koordinat& test_noktasi, const Koordinat& segment_ucu_1, const Koordinat& segment_ucu_2) {

	double min_x = min(segment_ucu_1.k_x, segment_ucu_2.k_x);
	double max_x = max(segment_ucu_1.k_x, segment_ucu_2.k_x);
	double min_y = min(segment_ucu_1.k_y, segment_ucu_2.k_y);
	double max_y = max(segment_ucu_1.k_y, segment_ucu_2.k_y);

	double hata_payi = 1e-5;

	return (test_noktasi.k_x >= min_x - hata_payi &&
		test_noktasi.k_x <= max_x + hata_payi &&
		test_noktasi.k_y >= min_y - hata_payi &&
		test_noktasi.k_y <= max_y + hata_payi);
}

bool SegmentlerKesisiyorMu(const Koordinat& p1, const Koordinat& p2,
	const Koordinat& p3, const Koordinat& p4,
	const DogruModeli& model_1, const DogruModeli& model_2)
{
	Koordinat kesisim;

	if (ModelKesisimBul(model_1, model_2, kesisim)) {

		if (NoktaSegmentUzerindeMi(kesisim, p1, p2) &&
			NoktaSegmentUzerindeMi(kesisim, p3, p4)) {
			return true;
		}
	}
	return false; // Paraleller veya kesişim segmentlerin dışında
}

void SiniflandirSegmentler(vector<SegmentBilgisi>& tum_segmentler)
{
	const Koordinat robot_orijin = { 0.0, 0.0 };

	for (size_t i = 0; i < tum_segmentler.size(); ++i) {
		auto& ana_segment = tum_segmentler[i];

		bool en_az_bir_nokta_gorunur = false;

		for (const auto& test_noktasi : ana_segment.ait_noktalar) {

			DogruModeli ray_model = NoktalardanModelOlustur(robot_orijin, test_noktasi);
			bool bu_nokta_engellendi = false;

			for (size_t j = 0; j < tum_segmentler.size(); ++j) {
				if (i == j) continue;

				const auto& engelleyici_segment = tum_segmentler[j];

				if (SegmentlerKesisiyorMu(robot_orijin, test_noktasi,
					engelleyici_segment.uc_1, engelleyici_segment.uc_2,
					ray_model, engelleyici_segment.model))
				{
					bu_nokta_engellendi = true;
					break;
				}
			}

			if (!bu_nokta_engellendi) {
				en_az_bir_nokta_gorunur = true;
				break;
			}
		}

		if (en_az_bir_nokta_gorunur) {
			ana_segment.gorunur = true;
		}
		else {
			ana_segment.gorunur = false;
		}
	}
}

vector<tuple<Koordinat, int, int>> KesisimAnalizcisi(const vector<SegmentBilgisi>& segment_listesi) {

	vector<tuple<Koordinat, int, int>> gecerli_kesisimler;

	for (size_t i = 0; i < segment_listesi.size(); ++i) {
		for (size_t j = i + 1; j < segment_listesi.size(); ++j) {

			Koordinat kesisme_noktasi;

			if (ModelKesisimBul(segment_listesi[i].model, segment_listesi[j].model, kesisme_noktasi)) {

				const auto& segment_1 = segment_listesi[i];
				const auto& segment_2 = segment_listesi[j];

				if (NoktaSegmentUzerindeMi(kesisme_noktasi, segment_1.uc_1, segment_1.uc_2) &&
					NoktaSegmentUzerindeMi(kesisme_noktasi, segment_2.uc_1, segment_2.uc_2)) {

					gecerli_kesisimler.emplace_back(kesisme_noktasi, static_cast<int>(i), static_cast<int>(j));
				}
			}
		}
	}
	return gecerli_kesisimler;
}

#ifndef PI
#define PI 3.141592653589793
#endif

double ModellerArasiAciBul(const DogruModeli& model_1, const DogruModeli& model_2) {
	double egim_1 = model_1.kat_A;
	double egim_2 = model_2.kat_A;
	double aci_radyan;

	if (isinf(egim_1) && isinf(egim_2)) {
		return 0.0;
	}
	else if (isinf(egim_1)) {
		aci_radyan = atan(fabs(1.0 / egim_2));
	}
	else if (isinf(egim_2)) {
		aci_radyan = atan(fabs(1.0 / egim_1));
	}
	else {
		aci_radyan = atan(fabs((egim_1 - egim_2) / (1.0 + egim_1 * egim_2)));
	}

	return aci_radyan * 180.0 / PI;
}

static ID3D11Device* g_pd3dDevice = nullptr;
static ID3D11DeviceContext* g_pd3dDeviceContext = nullptr;
static IDXGISwapChain* g_pSwapChain = nullptr;
static ID3D11RenderTargetView* g_mainRenderTargetView = nullptr;
bool CreateDeviceD3D(HWND hWnd);
void CleanupDeviceD3D();
void CreateRenderTarget();
void CleanupRenderTarget();
LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

// --- Ana Giriş Noktası ---
int WINAPI wWinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE hPrevInstance, _In_ LPWSTR lpCmdLine, _In_ int nCmdShow)
{
	// --- Veri Depolama Alanları ---
	vector<Koordinat> nokta_bulutu_verisi;
	vector<SegmentBilgisi> bulunan_segment_verisi;
	vector<tuple<Koordinat, int, int>> kesisen_segment_verisi;
	vector<double> kesim_mesafeleri;
	vector<double> kesim_acilari;

	vector<float> grafik_tum_x;
	vector<float> grafik_tum_y;

	// --- Veri İndirme ---
	const char* config_url = "https://gist.githubusercontent.com/NurayKBT/047778a7e4860a5a7aacec5e7616296b/raw/fa9d8acc8f30cb26371b3a75caaca59318bd49e0/ayarlar.toml";
	const char* yerel_config_adi = "config.toml";

	if (!YapilandirmaDosyasiniIndir(config_url, yerel_config_adi)) {
		return 1;
	}

	// --- Dosya Açma ve Ayrıştırma ---
	ifstream config_dosyasi(yerel_config_adi);
	if (!config_dosyasi.is_open()) {
		string hata_metni = "Indirilen dosya acilamadi: ";
		hata_metni += yerel_config_adi;
		MessageBoxA(NULL, hata_metni.c_str(), "Dosya Hatasi", MB_OK | MB_ICONERROR);
		return 1;
	}

	float param_aci_min = 0.0f;
	float param_aci_max = 0.0f;
	float param_aci_artis = 0.0f;
	float param_mesafe_min = 0.0f;
	float param_mesafe_max = 0.0f;
	vector<float> veri_mesafeler;
	vector<float> veri_yogunluklar;
	string okunan_satir;
	bool liste_okuma_aktif = false;
	vector<float>* gecerli_liste = nullptr;

	while (getline(config_dosyasi, okunan_satir)) {
		okunan_satir.erase(0, okunan_satir.find_first_not_of(" \t"));
		if (okunan_satir.empty() || okunan_satir[0] == '[') continue;

		if (okunan_satir.find("ranges") != string::npos) {
			liste_okuma_aktif = true;
			gecerli_liste = &veri_mesafeler;
			size_t parantezPos = okunan_satir.find('[');
			if (parantezPos != string::npos) okunan_satir = okunan_satir.substr(parantezPos + 1);
			else continue;
		}
		else if (okunan_satir.find("intensities") != string::npos) {
			liste_okuma_aktif = true;
			gecerli_liste = &veri_yogunluklar;
			size_t parantezPos = okunan_satir.find('[');
			if (parantezPos != string::npos) okunan_satir = okunan_satir.substr(parantezPos + 1);
			else continue;
		}

		if (liste_okuma_aktif) {
			if (okunan_satir.find(']') != string::npos) {
				liste_okuma_aktif = false;
				okunan_satir = okunan_satir.substr(0, okunan_satir.find(']'));
			}
			stringstream ss(okunan_satir);
			string okunan_deger;
			while (getline(ss, okunan_deger, ',')) {
				okunan_deger.erase(0, okunan_deger.find_first_not_of(" \t"));
				okunan_deger.erase(okunan_deger.find_last_not_of(" \t") + 1);
				if (!okunan_deger.empty() && gecerli_liste != nullptr) {
					try {
						gecerli_liste->push_back(stof(okunan_deger));
					}
					catch (...) {}
				}
			}
			if (!liste_okuma_aktif) continue;
		}

		size_t esitPos = okunan_satir.find('=');
		if (esitPos == string::npos) continue;
		string anahtar = okunan_satir.substr(0, esitPos);
		string deger = okunan_satir.substr(esitPos + 1);
		anahtar.erase(0, anahtar.find_first_not_of(" \t"));
		anahtar.erase(anahtar.find_last_not_of(" \t") + 1);
		deger.erase(0, deger.find_first_not_of(" \t"));
		deger.erase(deger.find_last_not_of(" \t") + 1);
		size_t yorumPos = deger.find('#');
		if (yorumPos != string::npos)
			deger = deger.substr(0, yorumPos);
		if (!deger.empty() && deger.front() == '"')
			deger.erase(0, 1);
		if (!deger.empty() && deger.back() == '"')
			deger.pop_back();

		try {
			if (anahtar == "angle_min") param_aci_min = stof(deger);
			else if (anahtar == "angle_max") param_aci_max = stof(deger);
			else if (anahtar == "angle_increment") param_aci_artis = stof(deger);
			else if (anahtar == "range_min") param_mesafe_min = stof(deger);
			else if (anahtar == "range_max") param_mesafe_max = stof(deger);
		}
		catch (...) {}
	}
	config_dosyasi.close();

	// --- Hesaplama Aşamaları ---
	nokta_bulutu_verisi = KutupsalDonusum(veri_mesafeler, param_aci_min, param_aci_artis, param_mesafe_min, param_mesafe_max);
	bulunan_segment_verisi = SegmentleriBoslugaGoreBol(nokta_bulutu_verisi);
	SiniflandirSegmentler(bulunan_segment_verisi);

	vector<SegmentBilgisi> gorunur_segmentler;
	vector<SegmentBilgisi> engellenmis_segmentler;
	for (const auto& segment : bulunan_segment_verisi) {
		if (segment.gorunur) {
			gorunur_segmentler.push_back(segment);
		}
		else {
			engellenmis_segmentler.push_back(segment);
		}
	}

	int engellenen_sayisi = static_cast<int>(engellenmis_segmentler.size());

	auto gecici_kesisimler = KesisimAnalizcisi(gorunur_segmentler);

	for (const auto& [nokta, i, j] : gecici_kesisimler) {

		double aci = ModellerArasiAciBul(gorunur_segmentler[i].model, gorunur_segmentler[j].model);

		if (aci >= 45.0) {
			double robota_uzaklik = sqrt(pow(nokta.k_x, 2) + pow(nokta.k_y, 2));
			kesisen_segment_verisi.push_back({ nokta, i, j });
			kesim_mesafeleri.push_back(robota_uzaklik);
			kesim_acilari.push_back(aci);
		}
	}

	// --- Grafik için Veri Hazırlama ---
	for (const auto& nokta : nokta_bulutu_verisi) {
		grafik_tum_x.push_back(static_cast<float>(nokta.k_x));
		grafik_tum_y.push_back(static_cast<float>(nokta.k_y));
	}

	// --- PENCERE OLUŞTURMA VE IMGUI BAŞLATMA ---
	WNDCLASSEX wc = { sizeof(WNDCLASSEX), CS_CLASSDC, WndProc, 0L, 0L, GetModuleHandle(nullptr), nullptr, nullptr, nullptr, nullptr, _T("Lidar Projesi"), nullptr };
	::RegisterClassEx(&wc);
	HWND hwnd = ::CreateWindow(wc.lpszClassName, _T("BLM210 Proje 1 - Lidar Veri Gorsellestirme"), WS_OVERLAPPEDWINDOW, 100, 100, 1280, 800, nullptr, nullptr, wc.hInstance, nullptr);

	if (!CreateDeviceD3D(hwnd)) {
		CleanupDeviceD3D();
		::UnregisterClass(wc.lpszClassName, wc.hInstance);
		return 1;
	}

	::ShowWindow(hwnd, SW_SHOWDEFAULT);
	::UpdateWindow(hwnd);

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImPlot::CreateContext();

	ImGui::StyleColorsDark();

	ImGui_ImplWin32_Init(hwnd);
	ImGui_ImplDX11_Init(g_pd3dDevice, g_pd3dDeviceContext);

	bool dongu_sonlansin = false;


	// --- ANA UYGULAMA DÖNGÜSÜ ---
	while (!dongu_sonlansin)
	{
		MSG mesaj;
		while (::PeekMessage(&mesaj, nullptr, 0U, 0U, PM_REMOVE))
		{
			::TranslateMessage(&mesaj);
			::DispatchMessage(&mesaj);
			if (mesaj.message == WM_QUIT)
				dongu_sonlansin = true;
		}
		if (dongu_sonlansin)
			break;

		ImGui_ImplDX11_NewFrame();
		ImGui_ImplWin32_NewFrame();
		ImGui::NewFrame();

		ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f));
		ImGui::SetNextWindowSize(ImGui::GetIO().DisplaySize);
		ImGui::Begin("GrafikPenceresi", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse);

		if (ImPlot::BeginPlot("LIDAR Grafigi", ImVec2(-1, -1), ImPlotFlags_Equal)) {

			// 1. Katman: Tüm Noktalar
			ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 3.0f, ImVec4(0.6f, 0.8f, 1.0f, 0.0f), 1.0f, ImVec4(0.6f, 0.8f, 1.0f, 1.0f));
			ImPlot::PlotScatter("Gecerli Noktalar", grafik_tum_x.data(), grafik_tum_y.data(), static_cast<int>(grafik_tum_x.size()));

			// 2. Katman: Robot
			float robot_x[] = { 0.0f };
			float robot_y[] = { 0.0f };
			ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 8.0f, ImVec4(1.0f, 0.0f, 0.0f, 1.0f), IMPLOT_AUTO, ImVec4(1.0f, 0.0f, 0.0f, 1.0f));
			ImPlot::PlotScatter("Robot", robot_x, robot_y, 1);

			// 3. Katman: Görünür Segmentler
			ImVec4 parlak_yesil = ImVec4(0.1f, 0.8f, 0.1f, 1.0f);

			for (size_t i = 0; i < gorunur_segmentler.size(); ++i)
			{
				const auto& segment = gorunur_segmentler[i];

				float cizgi_x[] = { (float)segment.uc_1.k_x, (float)segment.uc_2.k_x };
				float cizgi_y[] = { (float)segment.uc_1.k_y, (float)segment.uc_2.k_y };

				vector<float> segment_nokta_x, segment_nokta_y;
				for (const auto& p : segment.ait_noktalar) {
					segment_nokta_x.push_back(static_cast<float>(p.k_x));
					segment_nokta_y.push_back(static_cast<float>(p.k_y));
				}

				char ss_cizgi_adi[128];
				char ss_nokta_adi[128];
				snprintf(ss_cizgi_adi, sizeof(ss_cizgi_adi), "d%zu (%.2fm)", i + 1, segment.mesafe);
				snprintf(ss_nokta_adi, sizeof(ss_nokta_adi), "d%zu noktalari (%zu nokta)", i + 1, segment.ait_noktalar.size());

				ImPlot::PushStyleColor(ImPlotCol_Line, parlak_yesil);
				ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 3.0f);
				ImPlot::PlotLine(ss_cizgi_adi, cizgi_x, cizgi_y, 2);
				ImPlot::PopStyleVar();
				ImPlot::PopStyleColor();

				ImPlot::PushStyleColor(ImPlotCol_MarkerFill, parlak_yesil);
				ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, parlak_yesil);
				ImPlot::PushStyleVar(ImPlotStyleVar_Marker, ImPlotMarker_Circle);
				ImPlot::PushStyleVar(ImPlotStyleVar_MarkerSize, 3.0f);
				ImPlot::PlotScatter(ss_nokta_adi, segment_nokta_x.data(), segment_nokta_y.data(), static_cast<int>(segment_nokta_x.size()));
				ImPlot::PopStyleVar(2);
				ImPlot::PopStyleColor(2);

				string segment_adi = "d" + to_string(i + 1);
				ImPlot::PlotText(segment_adi.c_str(), (segment.uc_1.k_x + segment.uc_2.k_x) / 2.0, (segment.uc_1.k_y + segment.uc_2.k_y) / 2.0);
			}

			// 3.1. KATMAN: Engellenen Segmentler
			ImVec4 soluk_gri_cizgi = ImVec4(0.5f, 0.5f, 0.5f, 0.3f);
			ImVec4 soluk_gri_nokta = ImVec4(0.5f, 0.5f, 0.5f, 0.3f);

			char cizgi_grup_adi[128];
			char nokta_grup_adi[128];
			snprintf(cizgi_grup_adi, sizeof(cizgi_grup_adi), "Engellenen Dogrular (%d)", engellenen_sayisi);
			snprintf(nokta_grup_adi, sizeof(nokta_grup_adi), "Engellenen Noktalar");

			for (size_t i = 0; i < engellenmis_segmentler.size(); ++i)
			{
				const auto& segment = engellenmis_segmentler[i];

				float cizgi_x[] = { (float)segment.uc_1.k_x, (float)segment.uc_2.k_x };
				float cizgi_y[] = { (float)segment.uc_1.k_y, (float)segment.uc_2.k_y };

				vector<float> segment_nokta_x, segment_nokta_y;
				for (const auto& p : segment.ait_noktalar) {
					segment_nokta_x.push_back(static_cast<float>(p.k_x));
					segment_nokta_y.push_back(static_cast<float>(p.k_y));
				}

				ImPlot::PushStyleColor(ImPlotCol_Line, soluk_gri_cizgi);
				ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 1.5f);
				ImPlot::PlotLine(cizgi_grup_adi, cizgi_x, cizgi_y, 2);
				ImPlot::PopStyleVar();
				ImPlot::PopStyleColor();

				ImPlot::PushStyleColor(ImPlotCol_MarkerFill, soluk_gri_nokta);
				ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, soluk_gri_nokta);
				ImPlot::PushStyleVar(ImPlotStyleVar_Marker, ImPlotMarker_Circle);
				ImPlot::PushStyleVar(ImPlotStyleVar_MarkerSize, 2.0f);
				ImPlot::PlotScatter(nokta_grup_adi, segment_nokta_x.data(), segment_nokta_y.data(), static_cast<int>(segment_nokta_x.size()));
				ImPlot::PopStyleVar(2);
				ImPlot::PopStyleColor(2);
			}

			// 4. Katman: Geçerli Kesişimler
			const char* kesisim_grup_adi = "Kesisim";
			const char* mesafe_grup_adi = "Mesafe Cizgisi";

			for (size_t k = 0; k < kesisen_segment_verisi.size(); ++k)
			{
				auto& [nokta, i, j] = kesisen_segment_verisi[k];
				double mesafe = kesim_mesafeleri[k];
				double aci = kesim_acilari[k];

				float kesisim_x[] = { (float)nokta.k_x };
				float kesisim_y[] = { (float)nokta.k_y };
				ImPlot::SetNextMarkerStyle(ImPlotMarker_Cross, 10.0f, ImVec4(1.0f, 1.0f, 0.0f, 1.0f), 3.0f, ImVec4(1.0f, 1.0f, 0.0f, 1.0f));
				ImPlot::PlotScatter(kesisim_grup_adi, kesisim_x, kesisim_y, 1);

				stringstream ss_aci;
				ss_aci << fixed << setprecision(0) << aci << "' (d" << (i + 1) << " n d" << (j + 1) << ")";
				ImGui::PushStyleColor(ImGuiCol_PopupBg, IM_COL32(255, 255, 255, 255));
				ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(0, 0, 0, 255));
				ImVec4 text_col = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);
				ImPlot::Annotation(nokta.k_x, nokta.k_y, text_col, ImVec2(10, -10), true, "%s", ss_aci.str().c_str());
				ImGui::PopStyleColor(2);

				float mesafe_cizgisi_x[] = { 0.0f, (float)nokta.k_x };
				float mesafe_cizgisi_y[] = { 0.0f, (float)nokta.k_y };
				ImPlot::SetNextLineStyle(ImVec4(1.0f, 0.0f, 0.0f, 1.0f), 2.0f);
				ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 2.0f);
				ImPlot::PlotLine(mesafe_grup_adi, mesafe_cizgisi_x, mesafe_cizgisi_y, 2);
				ImPlot::PopStyleVar();

				stringstream ss_mesafe;
				ss_mesafe << fixed << setprecision(2) << mesafe << "m";
				ImGui::PushStyleColor(ImGuiCol_PopupBg, IM_COL32(255, 255, 0, 255));
				ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(0, 0, 0, 255));
				ImPlot::Annotation(nokta.k_x / 2.0, nokta.k_y / 2.0, text_col, ImVec2(0, 0), true, "%s", ss_mesafe.str().c_str());
				ImGui::PopStyleColor(2);
			}

			ImPlot::EndPlot();
		}
		ImGui::End();

		// --- EKRANA BASMA (RENDERING) ---
		ImGui::Render();
		const float arka_plan_rengi[4] = { 0.1f, 0.1f, 0.1f, 1.00f };
		g_pd3dDeviceContext->OMSetRenderTargets(1, &g_mainRenderTargetView, nullptr);
		g_pd3dDeviceContext->ClearRenderTargetView(g_mainRenderTargetView, arka_plan_rengi);
		ImGui_ImplDX11_RenderDrawData(ImGui::GetDrawData());
		g_pSwapChain->Present(1, 0);
	}

	// --- TEMİZLİK ---
	ImGui_ImplDX11_Shutdown();
	ImGui_ImplWin32_Shutdown();
	ImPlot::DestroyContext();
	ImGui::DestroyContext();

	CleanupDeviceD3D();
	::DestroyWindow(hwnd);
	::UnregisterClass(wc.lpszClassName, wc.hInstance);

	return 0;
}

bool CreateDeviceD3D(HWND hWnd)
{
	DXGI_SWAP_CHAIN_DESC sd;
	ZeroMemory(&sd, sizeof(sd));
	sd.BufferCount = 2;
	sd.BufferDesc.Width = 0;
	sd.BufferDesc.Height = 0;
	sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
	sd.BufferDesc.RefreshRate.Numerator = 60;
	sd.BufferDesc.RefreshRate.Denominator = 1;
	sd.Flags = DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH;
	sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
	sd.OutputWindow = hWnd;
	sd.SampleDesc.Count = 1;
	sd.SampleDesc.Quality = 0;
	sd.Windowed = TRUE;
	sd.SwapEffect = DXGI_SWAP_EFFECT_DISCARD;

	UINT createDeviceFlags = 0;
	D3D_FEATURE_LEVEL featureLevel;
	const D3D_FEATURE_LEVEL featureLevelArray[2] = { D3D_FEATURE_LEVEL_11_0, D3D_FEATURE_LEVEL_10_0, };
	if (D3D11CreateDeviceAndSwapChain(nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, createDeviceFlags, featureLevelArray, 2, D3D11_SDK_VERSION, &sd, &g_pSwapChain, &g_pd3dDevice, &featureLevel, &g_pd3dDeviceContext) != S_OK)
		return false;

	CreateRenderTarget();
	return true;
}

void CleanupDeviceD3D()
{
	CleanupRenderTarget();
	if (g_pSwapChain) { g_pSwapChain->Release(); g_pSwapChain = nullptr; }
	if (g_pd3dDeviceContext) { g_pd3dDeviceContext->Release(); g_pd3dDeviceContext = nullptr; }
	if (g_pd3dDevice) { g_pd3dDevice->Release(); g_pd3dDevice = nullptr; }
}

void CreateRenderTarget()
{
	ID3D11Texture2D* pBackBuffer = nullptr;
	HRESULT hr = g_pSwapChain->GetBuffer(0, IID_PPV_ARGS(&pBackBuffer));
	if (SUCCEEDED(hr) && pBackBuffer != nullptr)
	{
		g_pd3dDevice->CreateRenderTargetView(pBackBuffer, nullptr, &g_mainRenderTargetView);
		pBackBuffer->Release();
	}
}

void CleanupRenderTarget()
{
	if (g_mainRenderTargetView) { g_mainRenderTargetView->Release(); g_mainRenderTargetView = nullptr; }
}

extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	if (ImGui_ImplWin32_WndProcHandler(hWnd, msg, wParam, lParam))
		return true;

	switch (msg)
	{
	case WM_SIZE:
		if (g_pd3dDevice != nullptr && wParam != SIZE_MINIMIZED)
		{
			CleanupRenderTarget();
			g_pSwapChain->ResizeBuffers(0, (UINT)LOWORD(lParam), (UINT)HIWORD(lParam), DXGI_FORMAT_UNKNOWN, 0);
			CreateRenderTarget();
		}
		return 0;
	case WM_SYSCOMMAND:
		if ((wParam & 0xfff0) == SC_KEYMENU)
			return 0;
		break;
	case WM_DESTROY:
		::PostQuitMessage(0);
		return 0;
	}
	return ::DefWindowProc(hWnd, msg, wParam, lParam);
}