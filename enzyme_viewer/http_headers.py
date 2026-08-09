from __future__ import annotations

import gzip
from typing import Iterable


DEFAULT_GZIP_MIMETYPES = {
    "application/json",
    "application/javascript",
    "application/xml",
    "image/svg+xml",
    "text/css",
    "text/html",
    "text/javascript",
    "text/plain",
    "text/xml",
}

def build_content_security_policy(script_nonce: str | None = None) -> str:
    script_src = "script-src 'self'"
    if script_nonce:
        script_src += f" 'nonce-{script_nonce}'"
    return (
        "default-src 'self'; "
        f"{script_src}; "
        "style-src 'self' 'unsafe-inline'; "
        "img-src 'self' data: blob:; "
        "font-src 'self' data:; "
        "connect-src 'self'; "
        "object-src 'none'; "
        "base-uri 'self'; "
        "frame-ancestors 'none'"
    )


DEFAULT_CONTENT_SECURITY_POLICY = build_content_security_policy()


def append_vary_header(response, value: str) -> None:
    current = response.headers.get("Vary", "")
    values = [item.strip() for item in current.split(",") if item.strip()]
    if value.lower() not in {item.lower() for item in values}:
        values.append(value)
        response.headers["Vary"] = ", ".join(values)


def gzip_response_if_appropriate(
    response,
    *,
    request,
    min_bytes: int,
    mimetypes: Iterable[str] = DEFAULT_GZIP_MIMETYPES,
):
    accept_encoding = request.headers.get("Accept-Encoding", "").lower()
    if "gzip" not in accept_encoding:
        return response
    if request.method == "HEAD":
        return response
    if response.status_code < 200 or response.status_code >= 300:
        return response
    if response.status_code in {204, 304}:
        return response
    if response.headers.get("Content-Encoding"):
        return response
    if response.headers.get("ETag"):
        return response
    if response.direct_passthrough or response.is_streamed:
        return response
    if response.mimetype not in set(mimetypes):
        return response

    body = response.get_data()
    if len(body) < min_bytes:
        return response

    compressed = gzip.compress(body)
    if len(compressed) >= len(body):
        return response

    response.set_data(compressed)
    response.headers["Content-Encoding"] = "gzip"
    append_vary_header(response, "Accept-Encoding")
    return response


def apply_standard_response_headers(
    response,
    *,
    request,
    gzip_min_bytes: int,
    script_nonce: str | None = None,
):
    response.headers.setdefault(
        "Content-Security-Policy",
        build_content_security_policy(script_nonce),
    )
    response.headers.setdefault("X-Content-Type-Options", "nosniff")
    response.headers.setdefault("Referrer-Policy", "same-origin")
    response.headers.setdefault(
        "Permissions-Policy",
        "geolocation=(), microphone=(), camera=()",
    )
    return gzip_response_if_appropriate(
        response,
        request=request,
        min_bytes=gzip_min_bytes,
    )
