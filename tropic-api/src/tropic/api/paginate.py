"""Pagination configuration for FastAPI."""

from typing import TypeVar

from fastapi import Query
from fastapi_pagination import Page
from fastapi_pagination.customization import CustomizedPage, UseParamsFields

T = TypeVar("T")

BigPage = CustomizedPage[
    Page[T],
    UseParamsFields(
        size=Query(100, ge=1, le=1000),
    ),
]
