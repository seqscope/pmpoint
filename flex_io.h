#ifndef __FLEXIO_H__
#define __FLEXIO_H__

#include <string>
#include <cstdint>
#include <memory>
#include <cstdio>
#include <curl/curl.h>
#include <thread>
#include <chrono>

class FlexReader {
public:
    virtual ~FlexReader() = default;
    virtual bool open(const char* filename) = 0;
    //virtual bool read(uint64_t offset, uint64_t length, std::string& buffer) = 0;
    virtual bool read_at(uint64_t offset, uint64_t length, std::string& buffer) = 0;
    virtual uint64_t size_hint() const { return 0; }
    virtual bool is_open() const = 0;    
    virtual void close() = 0;
};

class FlexFileReader : public FlexReader {
    std::string path_;
    FILE* fp_ = nullptr;
    uint64_t size_ = 0;

public:
    FlexFileReader();
    ~FlexFileReader() override;
    bool open(const char* uri) override;
    bool read_at(uint64_t offset, uint64_t length, std::string& buffer) override;
    uint64_t size_hint() const override { return size_; }
    bool is_open() const override { return fp_ != nullptr; }
    void close() override { if ( is_open() ) { fclose(fp_); fp_ = nullptr; } }
};

class FlexHttpReader : public FlexReader {
    std::string url_;
    CURL* curl_ = nullptr;
    uint64_t size_ = 0;
    bool is_open_ = false;
    static size_t write_to_string(void* p, size_t sz, size_t nm, void* ud);
    bool parse_head();

public:
    FlexHttpReader();
    ~FlexHttpReader() override; // { if (curl_) { curl_easy_cleanup(curl_); curl_ = nullptr; } }
    bool open(const char* uri) override;
    bool read_at(uint64_t offset, uint64_t length, std::string& buffer) override;
    uint64_t size_hint() const override { return size_; }
    bool is_open() const override { return curl_ != nullptr; }
    void close() override { if ( is_open() ) { curl_easy_cleanup(curl_); curl_ = nullptr; } }
};

class FlexReaderFactory {
public:
    static std::unique_ptr<FlexReader> create_reader(const char* path);
};

#endif