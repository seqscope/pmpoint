#include "flex_io.h"
#include "qgenlib/qgen_error.h"

std::unique_ptr<FlexReader> FlexReaderFactory::create_reader(const char* path) {
    std::unique_ptr<FlexReader> reader;
    std::string path_str(path);

    if ( path_str.rfind("s3://", 0) == 0 ) {
        reader = std::make_unique<FlexHttpReader>();
        const std::string prefix = "s3://";
        size_t slash_pos = path_str.find('/', prefix.length());
        if (slash_pos == std::string::npos) return nullptr; // Invalid S3 URL
        
        std::string bucket = path_str.substr(prefix.length(), slash_pos - prefix.length());
        std::string key = path_str.substr(slash_pos + 1);
        path_str = "https://" + bucket + ".s3.amazonaws.com/" + key;
        //notice("path_str = %s", path_str.c_str());
    }
    else if (path_str.rfind("http://", 0) == 0 || path_str.rfind("https://", 0) == 0) {
        reader = std::make_unique<FlexHttpReader>();
    } else {
        reader = std::make_unique<FlexFileReader>();
    }

    if (reader && reader->open(path_str.c_str())) {
        return reader;
    }
    else {
        error("Failed to open reader for path: %s", path_str.c_str());
        return nullptr;
    }
}

FlexFileReader::FlexFileReader() {}
FlexFileReader::~FlexFileReader() { close(); }

bool FlexFileReader::open(const char* uri) {
    path_.assign(uri);
    fp_ = std::fopen(path_.c_str(), "rb");
    if ( !fp_ ) {
        error("Failed to open file %s", uri);
        return false;
    }
    std::fseek(fp_, 0, SEEK_END);
    size_ = std::ftell(fp_);
    std::fseek(fp_, 0, SEEK_SET);
    return true;
}

bool FlexFileReader::read_at(uint64_t offset, uint64_t length, std::string& buffer)  {
    if ( !is_open() ) return false;
    if (std::fseek(fp_, static_cast<long>(offset), SEEK_SET) != 0) return false;
    buffer.resize(length);
    size_t got = std::fread((void*)buffer.c_str(), 1, length, fp_);
    if ( got != length ) {
        buffer.clear();
        return false;
    }
    return true;
}

bool FlexHttpReader::parse_head() {
    CURL* c = curl_easy_init();
    if (!c) return false;
    curl_easy_setopt(c, CURLOPT_URL, url_.c_str());
    curl_easy_setopt(c, CURLOPT_NOBODY, 1L);
    curl_easy_setopt(c, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(c, CURLOPT_MAXREDIRS, 5L);
    auto rc = curl_easy_perform(c);
    long code = 0; curl_easy_getinfo(c, CURLINFO_RESPONSE_CODE, &code);
    double cl = 0; curl_easy_getinfo(c, CURLINFO_CONTENT_LENGTH_DOWNLOAD, &cl);
    curl_easy_cleanup(c);
    if (rc != CURLE_OK) return false;
    if (code == 200 || code == 204) { if (cl > 0) size_ = static_cast<uint64_t>(cl); return true; }
    return false;
}

size_t FlexHttpReader::write_to_string(void* p, size_t sz, size_t nm, void* ud) {
    //notice("write_to_string called with sz=%zu, nm=%zu", sz, nm);
    auto* s = static_cast<std::string*>(ud);
    s->append(static_cast<char*>(p), sz * nm);
    return sz * nm;

    // size_t new_length = size * nmemb;
    // try { s->append((char*)contents, new_length); } 
    // catch (std::bad_alloc& e) { return 0; }
    // return new_length;
}    

FlexHttpReader::FlexHttpReader() {
    curl_global_init(CURL_GLOBAL_DEFAULT);
    curl_ = curl_easy_init();
}

FlexHttpReader::~FlexHttpReader() {
    close();
    curl_global_cleanup();
}

bool FlexHttpReader::open(const char* uri) {
    //curl_ = curl_easy_init();
    // if (!curl_) return false;
    // curl_easy_setopt(curl_, CURLOPT_FOLLOWLOCATION, 1L);
    // curl_easy_setopt(curl_, CURLOPT_MAXREDIRS, 5L);
    // curl_easy_setopt(curl_, CURLOPT_TCP_KEEPALIVE, 1L);
    // parse_head(); // optional
    // is_open_ = true;
    // return true;

    if (!curl_) return false;
    url_.assign(uri);
    curl_easy_setopt(curl_, CURLOPT_URL, url_.c_str());
    curl_easy_setopt(curl_, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl_, CURLOPT_UNRESTRICTED_AUTH, 1L); 
    curl_easy_setopt(curl_, CURLOPT_NOBODY, 1L); // HEAD request
    if (curl_easy_perform(curl_) == CURLE_OK) {
        curl_off_t length = -1;
        curl_easy_getinfo(curl_, CURLINFO_CONTENT_LENGTH_DOWNLOAD_T, &length);
        if (length > 0) {
            size_ = static_cast<uint64_t>(length);
            is_open_ = true;
        }
    }
    curl_easy_setopt(curl_, CURLOPT_NOBODY, 0L); // Reset for future GETs
    return is_open_;
}

bool FlexHttpReader::read_at(uint64_t offset, uint64_t length, std::string& buffer) {
    if (!curl_) return false;
    //std::string range = "bytes=" + std::to_string(offset) + "-" + std::to_string(offset + length - 1);
    std::string range = std::to_string(offset) + "-" + std::to_string(offset + length - 1);
    //notice("Requesting url: %s", url_.c_str());
    //notice("Requesting range: %s", range.c_str());
    buffer.clear();
    curl_easy_setopt(curl_, CURLOPT_URL, url_.c_str());
    curl_easy_setopt(curl_, CURLOPT_RANGE, range.c_str());
    curl_easy_setopt(curl_, CURLOPT_WRITEFUNCTION, write_to_string);
    curl_easy_setopt(curl_, CURLOPT_WRITEDATA, &buffer);

    int32_t max_attempt = 4;
    for (int32_t attempt = 0; attempt < max_attempt; ++attempt) {
        auto rc = curl_easy_perform(curl_);
        long code = 0; 
        curl_easy_getinfo(curl_, CURLINFO_RESPONSE_CODE, &code);
        //notice("attempt %d: code = %ld", attempt, code);
        if (rc == CURLE_OK && (code == 206 || (code == 200 && length == buffer.size())) && buffer.size() == length) return true;
        //notice("buffer.size() = %zu, expected length = %zu", buffer.size(), length);
        if (code == 429 || code == 503) { std::this_thread::sleep_for(std::chrono::milliseconds(150 << attempt)); continue; }
        break;
    }
    buffer.clear();
    return false;
}