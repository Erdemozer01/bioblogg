from django.urls import path
from .views import *

app_name = "blog"

urlpatterns = [
    path('', BlogHomeView.as_view(), name="anasayfa"),
    path('categories/', CategoriesView.as_view(), name="categories"),
    path('category/<slug:slug>/', CategoryView.as_view(), name="category"),
    path('<slug:category>/<slug:slug>/<int:pk>/<slug:author>-<int:author_id>/<created>/',
         PostDetailView.as_view(), name="post_detail"),
    path('post-delete/<pk>/', PostDeleteView.as_view(), name="post_delete"),
    path('category-delete/<pk>/', CategoryDeleteView.as_view(), name="category_delete"),
    path('comment/<pk>/', CommentDetailView.as_view(), name="comment_detail"),
    path('delete/<pk>/', comment_delete, name="delete_comment"),
    path('like/<id>/', comment_like, name="comment_like"),
    path('dislike/<id>/', comment_dislike, name="comment_dislike"),
    path('report/<id>/', report_comment, name="report_comment"),
    path('okundu/<pk>/', comment_read, name="comment_read"),
    path('begendim/<int:pk>/<slug:slug>/', like_post, name="like_post"),
    path('begenmedim/<pk>/<title>/', dislike_post, name="dislike_post"),
    path('comment-update/<pk>/<commentator>/', CommentUpdate.as_view(), name="comment_update"),
    path('profile/<username>/<pk>/', profile_view, name="profile"),
    path('profile-update/<user>/<pk>/', ProfileUpdateViewNonStaff.as_view(), name="profile_update"),
    path('iletişim/', blog_contact, name="blog_contact"),
    path('hakkımda/', about, name="about"),
]
